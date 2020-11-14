#include "GridBasedFluid.h"

#include "Utilities/IO.h"
#include "Utilities/Yaml.h"

namespace PhysX {

template <int Dim>
GridBasedFluid<Dim>::GridBasedFluid(const StaggeredGrid<Dim> &grid) :
	_grid(grid),
	_velocity(&_grid)
{
	_velocity.parallelForEach([&](const int axis, const VectorDi &face) {
			const VectorDr pos = _velocity[axis].position(face);
			if constexpr (Dim == 2) _velocity[axis][face] = axis == 0 ? -pos.y() : pos.x();
		});

	_advection = std::move(std::make_unique<SemiLagrangianAdvection<Dim>>());
}

template<int Dim>
void GridBasedFluid<Dim>::writeDescription(std::ofstream &fout) const
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

template<int Dim>
void GridBasedFluid<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
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

template<int Dim>
void GridBasedFluid<Dim>::saveFrame(const std::string &frameDir) const
{
	std::ofstream fout(frameDir + "/velocity.sav", std::ios::binary);
	_velocity.write(fout);
}

template<int Dim>
void GridBasedFluid<Dim>::loadFrame(const std::string &frameDir)
{
	std::ifstream fin(frameDir + "/velocity.sav", std::ios::binary);
	_velocity.read(fin);
}

template <int Dim>
void GridBasedFluid<Dim>::advance(const real dt)
{
	_advection->advect(_velocity, _velocity, dt);
}

template class GridBasedFluid<2>;
template class GridBasedFluid<3>;

}
