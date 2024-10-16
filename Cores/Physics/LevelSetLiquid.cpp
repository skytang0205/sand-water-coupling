#include "LevelSetLiquid.h"

#include "Geometries/SurfaceMesh.h"
#include "Utilities/Constants.h"
#include "Utilities/IO.h"
#include "Utilities/MathFunc.h"
#include "Utilities/Yaml.h"

namespace PhysX {

template <int Dim>
LevelSetLiquid<Dim>::LevelSetLiquid(const StaggeredGrid<Dim> &grid) :
	EulerianFluid<Dim>(grid),
	_levelSet(_grid.cellGrid()),
	_levelSetReinitializer(std::make_unique<FastMarchingReinitializer<Dim>>(_grid.cellGrid()))
{ }

template <int Dim>
void LevelSetLiquid<Dim>::writeDescription(YAML::Node &root) const
{
	EulerianFluid<Dim>::writeDescription(root);
	{ // Description of liquid.
		YAML::Node node;
		node["name"] = "liquid";
		node["data_mode"] = "dynamic";

		if constexpr (Dim == 2) {
			node["primitive_type"] = "line_list";
			node["material"]["diffuse_albedo"] = Vector4f(0, 0, 1, 1);
		}
		else {
			node["primitive_type"] = "triangle_list";
			node["material"]["diffuse_albedo"] = (Vector4f(147, 213, 220, 255) / 255).eval(); // Qingshui Blue
		}

		node["indexed"] = true;
		root["objects"].push_back(node);
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	EulerianFluid<Dim>::writeFrame(frameDir, staticDraw);
	{ // Write liquid.
		std::ofstream fout(frameDir + "/liquid.mesh", std::ios::binary);
		SurfaceMesh<Dim> liquidMesh(_levelSet);
		IO::writeValue(fout, uint(liquidMesh.positions.size()));
		for (const auto &pos : liquidMesh.positions)
			IO::writeValue(fout, pos.template cast<float>().eval());
		if constexpr (Dim == 3) {
			for (const auto &normal : liquidMesh.normals)
				IO::writeValue(fout, normal.template cast<float>().eval());
		}
		IO::writeValue(fout, uint(liquidMesh.indices.size()));
		IO::writeArray(fout, liquidMesh.indices.data(), liquidMesh.indices.size());
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::saveFrame(const std::string &frameDir) const
{
	EulerianFluid<Dim>::saveFrame(frameDir);
	std::ofstream fout(frameDir + "/liquidSdf.sav", std::ios::binary);
	_levelSet.signedDistanceField().save(fout);
}

template <int Dim>
void LevelSetLiquid<Dim>::loadFrame(const std::string &frameDir)
{
	EulerianFluid<Dim>::loadFrame(frameDir);
	std::ifstream fin(frameDir + "/liquidSdf.sav", std::ios::binary);
	_levelSet.signedDistanceField().load(fin);
}

template <int Dim>
void LevelSetLiquid<Dim>::initialize()
{
	reinitializeLevelSet();
	EulerianFluid<Dim>::initialize();
}

template <int Dim>
void LevelSetLiquid<Dim>::advectFields(const real dt)
{
	_advector->advect(_levelSet.signedDistanceField(), _velocity, dt);
	reinitializeLevelSet();

	EulerianFluid<Dim>::advectFields(dt);
}

template <int Dim>
void LevelSetLiquid<Dim>::applyBodyForces(const real dt)
{
	if (_enableGravity) {
		_velocity[1].parallelForEach([&](const VectorDi &face) {
			_velocity[1][face] -= real(kGravity) * dt;
		});
	}

	EulerianFluid<Dim>::applyBodyForces(dt);
}

template <int Dim>
void LevelSetLiquid<Dim>::projectVelocity(const real dt)
{
	if (_enableSurfaceTension && dt)
		_projector->project(_velocity, _boundaryHelper->fraction(), _boundaryHelper->velocity(), _levelSet, _surfaceTensionCoefficient * dt / _density * _velocity.invSpacing());
	else
		_projector->project(_velocity, _boundaryHelper->fraction(), _boundaryHelper->velocity(), _levelSet);

	_boundaryHelper->extrapolate(_velocity, _levelSet, _kExtrapMaxSteps);
	_boundaryHelper->enforce(_velocity);
}

template <int Dim>
void LevelSetLiquid<Dim>::reinitializeLevelSet()
{
	_levelSetReinitializer->reinitialize(_levelSet, _kLsReinitMaxSteps);
}

template class LevelSetLiquid<2>;
template class LevelSetLiquid<3>;

}
