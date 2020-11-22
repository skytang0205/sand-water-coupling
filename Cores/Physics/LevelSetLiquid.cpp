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
	_projector->project(_velocity, _fluidFraction, _levelSet.signedDistanceField());
	extrapolateVelocity();
}

template <int Dim>
void LevelSetLiquid<Dim>::extrapolateVelocity()
{
	const auto &liquidSdf = _levelSet.signedDistanceField();
	const auto isLiquidFace = [&](const int axis, const VectorDi &face)->bool {
		const VectorDi cell0 = face - VectorDi::Unit(axis);
		const VectorDi cell1 = face;
		return (liquidSdf.isValid(cell0) && Surface<Dim>::isInside(liquidSdf[cell0]))
			|| (liquidSdf.isValid(cell1) && Surface<Dim>::isInside(liquidSdf[cell1]));
	};

	auto newVelocity = _velocity;
	newVelocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (isLiquidFace(axis, face)) return;
		int cnt = 0;
		real sum = 0;
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi &nbFace = Grid<Dim>::neighbor(face, i);
			if (newVelocity[axis].isValid(nbFace) && isLiquidFace(axis, nbFace))
				sum += newVelocity[axis][nbFace], cnt++;
		}
		if (cnt > 0) newVelocity[axis][face] = sum / cnt;
	});

	const real bandWidth = _kExtrapMaxIters * _grid.spacing();
	_velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		const VectorDr pos = _velocity[axis].position(face);
		if (_levelSet.signedDistance(pos) > 0) {
			_velocity[axis][face] =
				bandWidth < 0 || _levelSet.signedDistance(pos) < bandWidth ? newVelocity[axis](_levelSet.closestPosition(pos)) : real(0);
		}
	});

	enforceBoundaryConditions();
}

template <int Dim>
void LevelSetLiquid<Dim>::reinitializeLevelSet()
{
	_levelSetReinitializer->reinitialize(_levelSet, _kLsReinitMaxIters);
}

template class LevelSetLiquid<2>;
template class LevelSetLiquid<3>;

}
