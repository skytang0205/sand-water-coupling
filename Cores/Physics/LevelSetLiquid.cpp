#include "LevelSetLiquid.h"

#include "Utilities/Constants.h"
#include "Utilities/IO.h"
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
	{ // Description of liquidSdf.
		YAML::Node node;
		node["name"] = "liquidSdf";
		node["data_mode"] = "semi-dynamic";
		node["primitive_type"] = "triangle_list";
		node["indexed"] = true;
		node["color_map"]["enabled"] = true;
		root["objects"].push_back(node);
	}
}

template <int Dim>
void LevelSetLiquid<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	EulerianFluid<Dim>::writeFrame(frameDir, staticDraw);
	{ // Write liquidSdf.
		std::ofstream fout(frameDir + "/liquidSdf.mesh", std::ios::binary);
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
			const float liquidSdf = float(_levelSet.signedDistanceField()[cell]);
			IO::writeValue(fout, liquidSdf);
			IO::writeValue(fout, liquidSdf);
			IO::writeValue(fout, liquidSdf);
			IO::writeValue(fout, liquidSdf);
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
			_velocity[1][face] -= kGravity * dt;
		});
	}

	EulerianFluid<Dim>::applyBodyForces(dt);
}

template <int Dim>
void LevelSetLiquid<Dim>::projectVelocity()
{
	_boundary->enforce(_velocity);
	_projector->project(_velocity, _boundary->fraction(), _boundary->velocity(), _levelSet.signedDistanceField());
	_boundary->extrapolate(_velocity, _levelSet, _kExtrapMaxIters);
}

template <int Dim>
void LevelSetLiquid<Dim>::reinitializeLevelSet()
{
	_levelSetReinitializer->reinitialize(_levelSet, _kLsReinitMaxIters);
}

template class LevelSetLiquid<2>;
template class LevelSetLiquid<3>;

}
