#include "EulerianFluid.h"

#include "Utilities/IO.h"
#include "Utilities/Yaml.h"

namespace PhysX {

template <int Dim>
EulerianFluid<Dim>::EulerianFluid(const StaggeredGrid<Dim> &grid) :
	_grid(grid),
	_velocity(&_grid),
	_fluidFraction(&_grid),
	_advector(std::make_unique<SemiLagrangianAdvector<Dim>>()),
	_projector(std::make_unique<EulerianProjector<Dim>>(_grid.cellGrid()))
{
	_velocity.parallelForEach([&](const int axis, const VectorDi &face) {
			const VectorDr pos = _velocity[axis].position(face);
			if constexpr (Dim == 2) _velocity[axis][face] = axis == 0 ? -pos.y() : pos.x();
		});
}

template <int Dim>
void EulerianFluid<Dim>::writeDescription(std::ofstream &fout) const
{
	YAML::Node root;
	root["dimension"] = Dim;
	// Description of fluid.
	{
		YAML::Node node;
		node["name"] = "fluid";
		node["data_mode"] = "semi-dynamic";
		node["primitive_type"] = "triangle_list";
		node["indexed"] = true;
		node["color_map"]["enabled"] = true;
		root["objects"].push_back(node);
	}
	fout << root << std::endl;
}

template <int Dim>
void EulerianFluid<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	{ // Write fluid.
		std::ofstream fout(frameDir + "/fluid.mesh", std::ios::binary);
		IO::writeValue(fout, uint(4 * _grid.cellCount()));
		const VectorDr a = VectorDr::Unit(0) * _grid.spacing() / 2;
		const VectorDr b = VectorDr::Unit(1) * _grid.spacing() / 2;
		_grid.forEachCell([&](const VectorDi &cell) {
			const VectorDr pos = _grid.cellCenter(cell);
			IO::writeValue(fout, (pos - a - b).cast<float>().eval());
			IO::writeValue(fout, (pos + a - b).cast<float>().eval());
			IO::writeValue(fout, (pos - a + b).cast<float>().eval());
			IO::writeValue(fout, (pos + a + b).cast<float>().eval());
		});
		_grid.forEachCell([&](const VectorDi &cell) {
			const VectorDr pos = _grid.cellCenter(cell);
			const float vel = float(_velocity(pos).norm());
			IO::writeValue(fout, vel);
			IO::writeValue(fout, vel);
			IO::writeValue(fout, vel);
			IO::writeValue(fout, vel);
		});
		if (staticDraw) {
			static constexpr uint indices[] = { 1, 2, 0, 1, 3, 2 };
			IO::writeValue(fout, uint(6 * _grid.cellCount()));
			for (int i = 0; i < _grid.cellCount(); i++) {
				for (int j = 0; j < 6; j++)
					IO::writeValue(fout, uint(i * 4 + indices[j]));
			}
		}
	}
}

template <int Dim>
void EulerianFluid<Dim>::saveFrame(const std::string &frameDir) const
{
	std::ofstream fout(frameDir + "/velocity.sav", std::ios::binary);
	_velocity.write(fout);
}

template <int Dim>
void EulerianFluid<Dim>::loadFrame(const std::string &frameDir)
{
	std::ifstream fin(frameDir + "/velocity.sav", std::ios::binary);
	_velocity.read(fin);
}

template <int Dim>
void EulerianFluid<Dim>::initialize()
{
	// TODO: check pipeline
	updateFluidFraction();
	enforceBoundaryConditions();
	projectVelocity();
}

template <int Dim>
void EulerianFluid<Dim>::advance(const real dt)
{
	// TODO: check pipeline
	updateColliders(dt);
	advectFields(dt);
	updateFluidFraction();
	projectVelocity();
}

template <int Dim>
void EulerianFluid<Dim>::advectFields(const real dt)
{
	_advector->advect(_velocity, _velocity, dt);
	enforceBoundaryConditions();
}

template <int Dim>
void EulerianFluid<Dim>::updateColliders(const real dt)
{
	bool dirty = false;
	for (const auto &collider : _colliders) {
		auto dynamicCollider = dynamic_cast<DynamicCollider<Dim> *>(collider.get());
		if (!dynamicCollider) {
			dirty = true;
		}
	}
	if (dirty) updateFluidFraction();
}

template <int Dim>
void EulerianFluid<Dim>::applyBodyForces(const real dt)
{
}

template <int Dim>
void EulerianFluid<Dim>::projectVelocity()
{
	_projector->project(_velocity, _fluidFraction);
	enforceBoundaryConditions();
}

template <int Dim>
void EulerianFluid<Dim>::updateFluidFraction()
{
	_fluidFraction.parallelForEach([&](const int axis, const VectorDi &face) {
		_fluidFraction[axis][face] = 1;
		if (!_fluidFraction.isBoundary(axis, face)) {
			for (const auto &collider : _colliders) {
				const VectorDi cell0 = face - VectorDi::Unit(axis);
				const VectorDi cell1 = face;
				_fluidFraction[axis][face] -= collider->surface()->fractionInside(_grid.cellCenter(cell0), _grid.cellCenter(cell1));
			}
		}
	});
}

template <int Dim>
void EulerianFluid<Dim>::extrapolateVelocity()
{
	static constexpr int kMaxIterations = 1;

	auto newVelocity = _velocity;
	auto visited = std::make_unique<StaggeredGridBasedData<Dim, uchar>>(&_grid);
	auto newVisited = std::make_unique<StaggeredGridBasedData<Dim, uchar>>(&_grid);
	visited->parallelForEach([&](const int axis, const VectorDi &face) {
		if (!((*visited)[axis][face] = _fluidFraction[axis][face] > 0))
			newVelocity[axis][face] = 0;
	});
	for (int iter = 0; iter < kMaxIterations; iter++) {
		newVelocity.parallelForEach([&](const int axis, const VectorDi &face) {
			if (!(*visited)[axis][face]) {
				int cnt = 0;
				real sum = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbFace = Grid<Dim>::neighbor(face, i);
					if (newVelocity[axis].isValid(nbFace) && (*visited)[axis][nbFace])
						sum += newVelocity[axis][nbFace], cnt++;
				}
				if (cnt > 0) {
					newVelocity[axis][face] = sum / cnt;
					(*newVisited)[axis][face] = 1;
				}
			}
			else (*newVisited)[axis][face] = 1;
		});
		visited.swap(newVisited);
	}
	_velocity = newVelocity;
}

template <int Dim>
void EulerianFluid<Dim>::enforceBoundaryConditions()
{
	extrapolateVelocity();

	_velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (_velocity.isBoundary(axis, face)) {
			int pos = int(_time * 50) % 100;;
			if (axis == 0 && face[0] == 0) _velocity[axis][face] = 10 * (face[1] == pos);
			if (axis == 0 && face[0] > 0) _velocity[axis][face] = real(10) / _grid.resolution().y();
			if (axis == 1) _velocity[axis][face] = 0;
		}
		for (const auto &collider : _colliders) {
			const VectorDr pos = _velocity[axis].position(face);
			if (collider->surface()->inside(pos)) {
				const VectorDr colliderVel = collider->velocityAt(pos);
				const VectorDr vel = _velocity(pos);
				const VectorDr n = collider->surface()->closestNormal(pos);
				if (n.any()) {
					const VectorDr velR = vel - colliderVel;
					VectorDr velT = velR - velR.dot(n) * n;
					if (const real mu = collider->frictionCoefficient(); mu && velT.any()) {
						const real velN = std::max(-velR.dot(n), real(0));
						velT *= std::max(1 - mu * velN / velT.norm(), real(0));
					}
					_velocity[axis][face] = (velT + colliderVel)[axis];
				}
				else _velocity[axis][face] = colliderVel[axis];
				break;
			}
		}
	});
}

template class EulerianFluid<2>;
template class EulerianFluid<3>;

}
