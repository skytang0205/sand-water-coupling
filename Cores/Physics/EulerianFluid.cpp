#include "EulerianFluid.h"

#include "Geometries/ImplicitSurface.h"
#include "Utilities/IO.h"
#include "Utilities/Yaml.h"

namespace PhysX {

template <int Dim>
EulerianFluid<Dim>::EulerianFluid(const StaggeredGrid<Dim> &grid) :
	_grid(grid),
	_velocity(&_grid),
	_advector(std::make_unique<SemiLagrangianAdvector<Dim>>()),
	_boundaryHelper(std::make_unique<EulerianBoundaryHelper<Dim>>(&_grid)),
	_projector(std::make_unique<EulerianProjector<Dim>>(_grid.cellGrid()))
{ }

template <int Dim>
void EulerianFluid<Dim>::writeDescription(YAML::Node &root) const
{
	{ // Description of neumann
		YAML::Node node;
		node["name"] = "neumann";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["indexed"] = false;
		node["material"]["diffuse_albedo"] = Vector4f(0.5f, 0.5f, 0.5, 1.0f);
		root["objects"].push_back(node);
	}
	if constexpr (Dim == 2) { // Description of velocity.
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
	{ // Write neumann.
		std::ofstream fout(frameDir + "/neumann.mesh", std::ios::binary);
		uint cnt = 0;
		const auto &boundaryFraction = _boundaryHelper->fraction();
		boundaryFraction.forEach([&](const int axis, const VectorDi &face) {
			if (!boundaryFraction.isBoundary(axis, face) && boundaryFraction[axis][face] == 1)
				cnt++;
		});
		IO::writeValue(fout, cnt);
		boundaryFraction.forEach([&](const int axis, const VectorDi &face) {
			if (!boundaryFraction.isBoundary(axis, face) && boundaryFraction[axis][face] == 1)
				IO::writeValue(fout, _grid.faceCenter(axis, face).template cast<float>().eval());
		});
		if constexpr (Dim == 3) {
			boundaryFraction.forEach([&](const int axis, const VectorDi &face) {
				if (!boundaryFraction.isBoundary(axis, face) && boundaryFraction[axis][face] == 1)
					IO::writeValue(fout, VectorDf::Unit(2).eval());
			});
		}
	}
	if constexpr (Dim == 2) { // Write velocity.
		std::ofstream fout(frameDir + "/velocity.mesh", std::ios::binary);
		IO::writeValue(fout, uint(2 * _grid.cellCount()));
		_grid.forEachCell([&](const VectorDi &cell) {
			const VectorDr pos = _grid.cellCenter(cell);
			const VectorDr dir = _velocity(pos).normalized() * _grid.spacing() * std::sqrt(real(Dim)) / 2;
			IO::writeValue(fout, pos.template cast<float>().eval());
			IO::writeValue(fout, (pos + dir).template cast<float>().eval());
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
	updateBoundary();
}

template <int Dim>
void EulerianFluid<Dim>::initialize()
{
	updateBoundary();
	projectVelocity();
}

template <int Dim>
void EulerianFluid<Dim>::advance(const real dt)
{
	updateColliders(dt);

	advectFields(dt);
	applyBodyForces(dt);
	projectVelocity(dt);
}

template <int Dim>
void EulerianFluid<Dim>::advectFields(const real dt)
{
	_advector->advect(_velocity, _velocity, dt);
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
	if (dirty) updateBoundary();
}

template <int Dim>
void EulerianFluid<Dim>::projectVelocity(const real dt)
{
	_projector->project(_velocity, _boundaryHelper->fraction(), _boundaryHelper->velocity());
	_boundaryHelper->extrapolate(_velocity, _kExtrapMaxSteps);
	_boundaryHelper->enforce(_velocity);
}

template <int Dim>
void EulerianFluid<Dim>::updateBoundary()
{
	_boundaryHelper->reset(_colliders, _domainBoundaryVelocity);
}

template class EulerianFluid<2>;
template class EulerianFluid<3>;

}
