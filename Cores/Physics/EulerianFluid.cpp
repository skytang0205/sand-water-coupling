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
	{ // Description of neumann
		YAML::Node node;
		node["name"] = "neumann";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["indexed"] = false;
		node["material"]["diffuse_albedo"] = Vector4f(0.5f, 0.5f, 0.5, 1.0f);
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
	{ // Write neumann.
		std::ofstream fout(frameDir + "/neumann.mesh", std::ios::binary);
		uint cnt = 0;
		_fluidFraction.forEach([&](const int axis, const VectorDi &face) {
			if (!_fluidFraction.isBoundary(axis, face) && _fluidFraction[axis][face] < 1) cnt++;
		});
		IO::writeValue(fout, cnt);
		_fluidFraction.forEach([&](const int axis, const VectorDi &face) {
			if (!_fluidFraction.isBoundary(axis, face) && _fluidFraction[axis][face] < 1)
				IO::writeValue(fout, _grid.faceCenter(axis, face).cast<float>().eval());
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
				const VectorDr pos0 = _grid.cellCenter(face - VectorDi::Unit(axis));
				const VectorDr pos1 = _grid.cellCenter(face);
				real phi0 = std::numeric_limits<real>::infinity();
				real phi1 = std::numeric_limits<real>::infinity();
				for (const auto &collider : _colliders) {
					phi0 = std::min(phi0, collider->surface()->signedDistance(pos0));
					phi1 = std::min(phi1, collider->surface()->signedDistance(pos1));
				}
				_fluidFraction[axis][face] -= Surface<Dim>::fraction(phi0, phi1);
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

	for (int iter = 0; iter < _kExtrapMaxIters || (_kExtrapMaxIters < 0 && visited->sum<size_t>() < visited->count()); iter++) {
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
		if (_domainBoundaryHandler && _velocity.isBoundary(axis, face))
			_velocity[axis][face] = _domainBoundaryHandler(axis, face);
		else
			for (const auto &collider : _colliders) {
				const VectorDr pos = _velocity[axis].position(face);
				if (collider->surface()->isInside(pos)) {
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
