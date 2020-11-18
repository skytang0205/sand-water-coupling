#include "LevelSetLiquid.h"

#include "Utilities/IO.h"
#include "Utilities/Yaml.h"

namespace PhysX {

template <int Dim>
LevelSetLiquid<Dim>::LevelSetLiquid(const StaggeredGrid<Dim> &grid) :
	EulerianFluid<Dim>(grid),
	_levelSet(_grid.cellGrid()),
	_levelSetReinitializer(std::make_unique<FastMarchingReinitializer<Dim>>(_grid.cellGrid(), -1))
{ }

template <int Dim>
void LevelSetLiquid<Dim>::writeDescription(std::ofstream &fout) const
{
	YAML::Node root;
	root["dimension"] = Dim;
	// Description of fluid.
	{
		YAML::Node node;
		node["name"] = "phi";
		node["data_mode"] = "semi-dynamic";
		node["primitive_type"] = "triangle_list";
		node["indexed"] = true;
		node["color_map"]["enabled"] = true;
		root["objects"].push_back(node);
	}
	fout << root << std::endl;
}

template <int Dim>
void LevelSetLiquid<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	{ // Write fluid.
		std::ofstream fout(frameDir + "/phi.mesh", std::ios::binary);
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
			const float phi = _levelSet.signedDistanceField()[cell];
			IO::writeValue(fout, phi);
			IO::writeValue(fout, phi);
			IO::writeValue(fout, phi);
			IO::writeValue(fout, phi);
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
	std::ofstream fout(frameDir + "/phi.sav", std::ios::binary);
	_levelSet.signedDistanceField().save(fout);
}

template <int Dim>
void LevelSetLiquid<Dim>::loadFrame(const std::string &frameDir)
{
	std::ifstream fin(frameDir + "/phi.sav", std::ios::binary);
	_levelSet.signedDistanceField().load(fin);
}

template <int Dim>
void LevelSetLiquid<Dim>::initialize()
{
}

template <int Dim>
void LevelSetLiquid<Dim>::advance(const real dt)
{
	_advector->advect(_levelSet.signedDistanceField(), _velocity, dt);
	_levelSetReinitializer->reinitialize(_levelSet);
}

template class LevelSetLiquid<2>;
template class LevelSetLiquid<3>;

}
