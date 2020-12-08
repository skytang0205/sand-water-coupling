#include "ParticleInCellLiquid.h"

#include <numbers>

namespace PhysX {

template <int Dim>
ParticleInCellLiquid<Dim>::ParticleInCellLiquid(const StaggeredGrid<Dim> &grid, const int particlesCntPerSubcell) :
	LevelSetLiquid<Dim>(grid),
	_particlesCntPerSubCell(particlesCntPerSubcell)
{ }

template <int Dim>
void ParticleInCellLiquid<Dim>::writeDescription(YAML::Node &root) const
{
	LevelSetLiquid<Dim>::writeDescription(root);
	if constexpr (Dim == 2) { // Description of particles.
		YAML::Node node;
		node["name"] = "particles";
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
	if constexpr (Dim == 2) { // Write particles.
		std::ofstream fout(frameDir + "/particles.mesh", std::ios::binary);
		IO::writeValue(fout, uint(_particles.size()));
		_particles.forEach([&](const int i) {
			IO::writeValue(fout, _particles.positions[i].template cast<float>().eval());
		});
	}
}

template <int Dim>
void ParticleInCellLiquid<Dim>::saveFrame(const std::string &frameDir) const
{
	EulerianFluid<Dim>::saveFrame(frameDir);
	{ // Save particles.
		std::ofstream fout(frameDir + "/particles.sav", std::ios::binary);
		IO::writeValue(fout, uint(_particles.size()));
		_particles.positions.save(fout);
	}
}

template <int Dim>
void ParticleInCellLiquid<Dim>::loadFrame(const std::string &frameDir)
{
	EulerianFluid<Dim>::loadFrame(frameDir);
	{ // Load particles.
		std::ifstream fin(frameDir + "/particles.sav", std::ios::binary);
		uint particlesCnt;
		IO::readValue(fin, particlesCnt);
		_particles.resize(particlesCnt);
		_particles.positions.load(fin);
		_particleVelocities.resize(&_particles);
		_particleVelocities.setZero();
	}
}

template <int Dim>
void ParticleInCellLiquid<Dim>::initialize()
{
	reinitializeParticles();
	LevelSetLiquid<Dim>::initialize();
}

template <int Dim>
void ParticleInCellLiquid<Dim>::advance(const real dt)
{
	updateColliders(dt);

	transferFromGridToParticles();
	advectFields(dt);
	applyParticleForces(dt);
	transferFromParticlesToGrid();

	applyBodyForces(dt);
	projectVelocity(dt);
}

template <int Dim>
void ParticleInCellLiquid<Dim>::advectFields(const real dt)
{
	_advector->advect(_particles, _velocity, dt);
	_boundaryHelper->enforce(_particles, _particleVelocities);
}

template <int Dim>
void ParticleInCellLiquid<Dim>::transferFromGridToParticles()
{
	_particles.parallelForEach([&](const int i) {
		const VectorDr pos = _particles.positions[i];
		_particleVelocities[i] = _velocity(pos);
	});
}

template <int Dim>
void ParticleInCellLiquid<Dim>::transferFromParticlesToGrid()
{
	StaggeredGridBasedScalarData<Dim> weightSum(_velocity.staggeredGrid());
	_velocity.setZero();

	transferFromParticlesToGrid(weightSum);

	_velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (weightSum[axis][face])
			_velocity[axis][face] /= weightSum[axis][face];
	});

	maintainGridBasedData(weightSum);
}

template <int Dim>
void ParticleInCellLiquid<Dim>::transferFromParticlesToGrid(StaggeredGridBasedScalarData<Dim> &weightSum)
{
	_particles.forEach([&](const int i) {
		const VectorDr pos = _particles.positions[i];
		const VectorDr vel = _particleVelocities[i];
		for (int axis = 0; axis < Dim; axis++) {
			for (const auto [face, weight] : _velocity[axis].grid()->linearIntrplDataPoints(pos)) {
				_velocity[axis][face] += vel[axis] * weight;
				weightSum[axis][face] += weight;
			}
		}
	});
}

template <int Dim>
void ParticleInCellLiquid<Dim>::maintainGridBasedData(StaggeredGridBasedScalarData<Dim> &weightSum)
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

	_particles.forEach([&](const int i) {
		const VectorDr pos = _particles.positions[i];
		const ImplicitSphere<Dim> sphere(pos, radius);
		for (const auto &cell : liquidSdf.grid()->cubicNearbyDataPoints(pos)) {
			if (liquidSdf.isValid(cell))
				liquidSdf[cell] = std::min(liquidSdf[cell], sphere.signedDistance(liquidSdf.position(cell)));
		}
	});

	_levelSetReinitializer->reinitialize(_levelSet, _kLsReinitMaxSteps);
}

template <int Dim>
void ParticleInCellLiquid<Dim>::reinitializeParticles()
{
	_particles.clear();

	auto &liquidSdf = _levelSet.signedDistanceField();
	const real dx = liquidSdf.spacing();
	const real radius = dx * real(1.1) / real(std::numbers::sqrt2);
	liquidSdf.forEach([&](const VectorDi &cell) {
		const VectorDr centerPos = liquidSdf.position(cell);
		for (int i = 0; i < (1 << Dim) * _particlesCntPerSubCell; i++) {
			const VectorDr pos = centerPos + VectorDr::Random() * dx * real(.5);
			if (_levelSet.signedDistance(pos) <= -radius)
				_particles.add(pos);
		}
	});
	_particleVelocities.resize(&_particles);
	_particleVelocities.setZero();
}

template class ParticleInCellLiquid<2>;
template class ParticleInCellLiquid<3>;

}
