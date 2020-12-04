#include "ParticleInCellLiquid.h"

#include <numbers>

namespace PhysX {

template <int Dim>
ParticleInCellLiquid<Dim>::ParticleInCellLiquid(const StaggeredGrid<Dim> &grid, const int markersCntPerSubcell) :
	LevelSetLiquid<Dim>(grid),
	_markersCntPerSubCell(markersCntPerSubcell)
{ }

template <int Dim>
void ParticleInCellLiquid<Dim>::writeDescription(YAML::Node &root) const
{
	LevelSetLiquid<Dim>::writeDescription(root);
	if constexpr (Dim == 2) { // Description of markers.
		YAML::Node node;
		node["name"] = "markers";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["material"]["diffuse_albedo"] = Vector4f(0, 0, 0, 1);
		node["indexed"] = false;
		root["objects"].push_back(node);
	}
}

template <int Dim>
void ParticleInCellLiquid<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	LevelSetLiquid<Dim>::writeFrame(frameDir, staticDraw);
	if constexpr (Dim == 2) { // Write Markers.
		std::ofstream fout(frameDir + "/markers.mesh", std::ios::binary);
		IO::writeValue(fout, uint(_markerPositions.size()));
		_markerPositions.forEach([&](const int i) {
			IO::writeValue(fout, _markerPositions[i].template cast<float>().eval());
		});
	}
}

template <int Dim>
void ParticleInCellLiquid<Dim>::saveFrame(const std::string &frameDir) const
{
	LevelSetLiquid<Dim>::saveFrame(frameDir);
	{ // Save markerPositions.
		std::ofstream fout(frameDir + "/markerPositions.sav", std::ios::binary);
		_markerPositions.save(fout);
	}
}

template <int Dim>
void ParticleInCellLiquid<Dim>::loadFrame(const std::string &frameDir)
{
	LevelSetLiquid<Dim>::loadFrame(frameDir);
	reinitializeMarkers();
	{ // Load markerPositions.
		std::ifstream fin(frameDir + "/markerPositions.sav", std::ios::binary);
		_markerPositions.load(fin);
	}
}

template <int Dim>
void ParticleInCellLiquid<Dim>::initialize()
{
	reinitializeMarkers();
	LevelSetLiquid<Dim>::initialize();
}

template <int Dim>
void ParticleInCellLiquid<Dim>::advance(const real dt)
{
	updateColliders(dt);

	transferFromGridsToParticles();
	advectFields(dt);
	applyMarkerForces(dt);
	transferFromParticlesToGrids();

	applyBodyForces(dt);
	projectVelocity(dt);
}

template <int Dim>
void ParticleInCellLiquid<Dim>::advectFields(const real dt)
{
	_advector->advect(_markerPositions, _velocity, dt);
	_boundaryHelper->enforce(_markerPositions, _markerVelocities);
}

template <int Dim>
void ParticleInCellLiquid<Dim>::transferFromGridsToParticles()
{
	_markerVelocities.parallelForEach([&](const int i) {
		const VectorDr pos = _markerPositions[i];
		_markerVelocities[i] = _velocity(pos);
	});
}

template <int Dim>
void ParticleInCellLiquid<Dim>::transferFromParticlesToGrids()
{
	StaggeredGridBasedData<Dim> weightSum(_velocity.staggeredGrid());
	_velocity.setZero();

	transferFromParticlesToGrids(weightSum);

	_velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (weightSum[axis][face])
			_velocity[axis][face] /= weightSum[axis][face];
	});

	maintainGridsBasedData(weightSum);
}

template <int Dim>
void ParticleInCellLiquid<Dim>::transferFromParticlesToGrids(StaggeredGridBasedData<Dim> &weightSum)
{
	_markerVelocities.forEach([&](const int i) {
		const VectorDr pos = _markerPositions[i];
		const VectorDr vel = _markerVelocities[i];
		for (int axis = 0; axis < Dim; axis++) {
			for (const auto [face, weight] : _velocity[axis].grid()->linearIntrplDataPoints(pos)) {
				_velocity[axis][face] += vel[axis] * weight;
				weightSum[axis][face] += weight;
			}
		}
	});
}

template <int Dim>
void ParticleInCellLiquid<Dim>::maintainGridsBasedData(StaggeredGridBasedData<Dim> &weightSum)
{
	reinitializeLevelSet();
	_boundaryHelper->extrapolate(_velocity, _levelSet, weightSum, _kExtrapMaxSteps);
}

template <int Dim>
void ParticleInCellLiquid<Dim>::reinitializeLevelSet()
{
	_levelSet.clear();

	auto &liquidSdf = _levelSet.signedDistanceField();
	const real radius = liquidSdf.spacing() * real(1.1) / real(std::numbers::sqrt2);

	_markerPositions.forEach([&](const int i) {
		const VectorDr pos = _markerPositions[i];
		const ImplicitSphere<Dim> sphere(pos, radius);
		for (const auto &cell : liquidSdf.grid()->cubicNearbyDataPoints(pos)) {
			if (liquidSdf.isValid(cell))
				liquidSdf[cell] = std::min(liquidSdf[cell], sphere.signedDistance(liquidSdf.position(cell)));
		}
	});

	_levelSetReinitializer->reinitialize(_levelSet, _kLsReinitMaxSteps);
}

template <int Dim>
void ParticleInCellLiquid<Dim>::reinitializeMarkers()
{
	_markerPositions.clear();

	auto &liquidSdf = _levelSet.signedDistanceField();
	const real dx = liquidSdf.spacing();
	const real radius = dx * real(1.1) / real(std::numbers::sqrt2);
	liquidSdf.forEach([&](const VectorDi &cell) {
		const VectorDr centerPos = liquidSdf.position(cell);
		for (int i = 0; i < (1 << Dim) * _markersCntPerSubCell; i++) {
			const VectorDr pos = centerPos + VectorDr::Random() * dx * real(0.5);
			if (_levelSet.signedDistance(pos) <= -radius)
				_markerPositions.add(pos);
		}
	});
	_markerVelocities.resize(_markerPositions.size());
	_markerVelocities.setZero();
}

template class ParticleInCellLiquid<2>;
template class ParticleInCellLiquid<3>;

}
