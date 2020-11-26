#include "EulerianFluid.h"

#include "Geometries/ImplicitSurface.h"
#include "Utilities/IO.h"
#include "Utilities/Yaml.h"

namespace PhysX {

template <int Dim>
EulerianFluid<Dim>::EulerianFluid(const StaggeredGrid<Dim> &grid) :
	_grid(grid),
	_velocity(&_grid),
	_boundary(std::make_unique<EulerianBoundary<Dim>>(&_grid)),
	_advector(std::make_unique<MacCormackAdvector<Dim>>()),
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
		if constexpr (Dim == 3) {
			_grid.forEachCell([&](const VectorDi &cell) {
				IO::writeValue(fout, VectorDf::Unit(2).eval());
				IO::writeValue(fout, VectorDf::Unit(2).eval());
			});
		}
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
		const auto &boundaryFraction = _boundary->fraction();
		boundaryFraction.forEach([&](const int axis, const VectorDi &face) {
			if (!boundaryFraction.isBoundary(axis, face) && boundaryFraction[axis][face] == 1)
				cnt++;
		});
		IO::writeValue(fout, cnt);
		boundaryFraction.forEach([&](const int axis, const VectorDi &face) {
			if (!boundaryFraction.isBoundary(axis, face) && boundaryFraction[axis][face] == 1)
				IO::writeValue(fout, _grid.faceCenter(axis, face).cast<float>().eval());
		});
		if constexpr (Dim == 3) {
			boundaryFraction.forEach([&](const int axis, const VectorDi &face) {
				if (!boundaryFraction.isBoundary(axis, face) && boundaryFraction[axis][face] == 1)
					IO::writeValue(fout, VectorDf::Unit(2).eval());
			});
		}
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
void EulerianFluid<Dim>::applyBodyForces(const real dt)
{
}

template <int Dim>
void EulerianFluid<Dim>::projectVelocity(const real dt)
{
	_boundary->enforce(_velocity);
	_projector->project(_velocity, _boundary->fraction(), _boundary->velocity());
	_boundary->extrapolate(_velocity, _kExtrapMaxSteps);
}

template <int Dim>
void EulerianFluid<Dim>::updateBoundary()
{
	_boundary->reset(_colliders, _domainBoundaryVelocity);
}

template class EulerianFluid<2>;
template class EulerianFluid<3>;

}
