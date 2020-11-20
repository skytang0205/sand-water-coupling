#include "EulerianFluid.h"

#include "Geometries/ImplicitSurface.h"
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
{ }

template <int Dim>
void EulerianFluid<Dim>::writeDescription(YAML::Node &root) const
{
	{ // Description of velocity.
		YAML::Node node;
		node["name"] = "velocity";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "line_list";
		node["indexed"] = false;
		node["color_map"]["enabled"] = true;
		root["objects"].push_back(node);
	}
}

template <int Dim>
void EulerianFluid<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	{ // Write velocity.
		std::ofstream fout(frameDir + "/velocity.mesh", std::ios::binary);
		IO::writeValue(fout, uint(2 * _grid.cellCount()));
		_grid.forEachCell([&](const VectorDi &cell) {
			const VectorDr pos = _grid.cellCenter(cell);
			const VectorDr dir = _velocity(pos).normalized() * _grid.spacing() * std::sqrt(real(Dim)) / 2;
			IO::writeValue(fout, pos.cast<float>().eval());
			IO::writeValue(fout, (pos + dir).cast<float>().eval());
		});
		_grid.forEachCell([&](const VectorDi &cell) {
			const VectorDr pos = _grid.cellCenter(cell);
			const float vel = float(_velocity(pos).norm());
			IO::writeValue(fout, vel);
			IO::writeValue(fout, vel);
		});
	}
}

template <int Dim>
void EulerianFluid<Dim>::saveFrame(const std::string &frameDir) const
{
	std::ofstream fout(frameDir + "/velocity.sav", std::ios::binary);
	_velocity.save(fout);
}

template <int Dim>
void EulerianFluid<Dim>::loadFrame(const std::string &frameDir)
{
	std::ifstream fin(frameDir + "/velocity.sav", std::ios::binary);
	_velocity.load(fin);
	updateFluidFraction();
}

template <int Dim>
void EulerianFluid<Dim>::initialize()
{
	updateFluidFraction();
	extrapolateVelocity();
	projectVelocity();
}

template <int Dim>
void EulerianFluid<Dim>::advance(const real dt)
{
	updateColliders(dt);

	advectFields(dt);
	applyBodyForces(dt);
	projectVelocity();
}

template <int Dim>
void EulerianFluid<Dim>::advectFields(const real dt)
{
	_advector->advect(_velocity, _velocity, dt);
	extrapolateVelocity();
}

template <int Dim>
void EulerianFluid<Dim>::updateColliders(const real dt)
{
	bool dirty = false;
	for (const auto &collider : _colliders) {
		auto dynamicCollider = dynamic_cast<DynamicCollider<Dim> *>(collider.get());
		if (dynamicCollider) {
			dirty = true;
		}
	}
	if (dirty) updateFluidFraction();
}

template <int Dim>
void EulerianFluid<Dim>::applyBodyForces(const real dt)
{
	extrapolateVelocity();
}

template <int Dim>
void EulerianFluid<Dim>::projectVelocity()
{
	_projector->project(_velocity, _fluidFraction);
	extrapolateVelocity();
}

template <int Dim>
void EulerianFluid<Dim>::updateFluidFraction()
{
	_fluidFraction.setConstant(1);
	if (!_colliders.empty()) {
		_fluidFraction.parallelForEach([&](const int axis, const VectorDi &face) {
			if (!_fluidFraction.isBoundary(axis, face)) {
				for (const auto &collider : _colliders) {
					const VectorDi cell0 = face - VectorDi::Unit(axis);
					const VectorDi cell1 = face;
					_fluidFraction[axis][face] -= Surface<Dim>::fraction(
						collider->surface()->signedDistance(_grid.cellCenter(cell0)),
						collider->surface()->signedDistance(_grid.cellCenter(cell1)));
				}
			}
		});
	}
}

template <int Dim>
void EulerianFluid<Dim>::extrapolateVelocity()
{
	auto newVelocity = _velocity;
	auto visited = std::make_unique<StaggeredGridBasedData<Dim, uchar>>(&_grid);
	auto newVisited = std::make_unique<StaggeredGridBasedData<Dim, uchar>>(&_grid);
	visited->parallelForEach([&](const int axis, const VectorDi &face) {
		if (!((*visited)[axis][face] = _fluidFraction[axis][face] > 0))
			newVelocity[axis][face] = 0;
	});
	for (int iter = 0; iter < _kExtrapMaxIters; iter++) {
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
					(*newVisited)[axis][face] = true;
				}
			}
			else (*newVisited)[axis][face] = true;
		});
		visited.swap(newVisited);
	}
	_velocity = newVelocity;

	enforceBoundaryConditions();
}

template <int Dim>
void EulerianFluid<Dim>::enforceBoundaryConditions()
{
	_velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (_velocity.isBoundary(axis, face) && _domainBoundaryHandler) {
			_velocity[axis][face] = _domainBoundaryHandler(axis, face);
		}
		else {
			for (const auto &collider : _colliders) {
				const VectorDr pos = _velocity[axis].position(face);
				if (collider->surface()->isInside(pos)) {
					_velocity[axis][face] = collider->velocityAt(pos)[axis];
					break;
				}
			}
		}
	});
}

template class EulerianFluid<2>;
template class EulerianFluid<3>;

}
