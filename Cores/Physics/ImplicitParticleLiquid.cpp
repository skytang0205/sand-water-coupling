#include "ImplicitParticleLiquid.h"

#include <numbers>

namespace PhysX {

template <int Dim>
ImplicitParticleLiquid<Dim>::ImplicitParticleLiquid(const StaggeredGrid<Dim> &grid, const int markersCntPerSubcell, const real propOfPic) :
	ParticleInCellLiquid<Dim>(grid, markersCntPerSubcell),
	_propOfPic(std::clamp(propOfPic, real(0), real(1))),
	_deltaVelocity(&_grid)
{ }

template <int Dim>
void ImplicitParticleLiquid<Dim>::saveFrame(const std::string &frameDir) const
{
	ParticleInCellLiquid<Dim>::saveFrame(frameDir);
	{ // Save particleVelocities.
		std::ofstream fout(frameDir + "/particleVelocities.sav", std::ios::binary);
		_particleVelocities.save(fout);
	}
	{ // Save delteVelocity.
		std::ofstream fout(frameDir + "/delteVelocity.sav", std::ios::binary);
		_deltaVelocity.save(fout);
	}
}

template <int Dim>
void ImplicitParticleLiquid<Dim>::loadFrame(const std::string &frameDir)
{
	ParticleInCellLiquid<Dim>::loadFrame(frameDir);
	{ // Load particleVelocities.
		std::ifstream fin(frameDir + "/particleVelocities.sav", std::ios::binary);
		_particleVelocities.load(fin);
	}
	{ // Load delteVelocity.
		std::ifstream fin(frameDir + "/delteVelocity.sav", std::ios::binary);
		_deltaVelocity.load(fin);
	}
}

template <int Dim>
void ImplicitParticleLiquid<Dim>::transferFromGridToParticles()
{
	_deltaVelocity.parallelForEach([&](const int axis, const VectorDi &face) {
		_deltaVelocity[axis][face] = _velocity[axis][face] - _deltaVelocity[axis][face];
	});
	_particles.parallelForEach([&](const int i) {
		const VectorDr pos = _particles.positions[i];
		_particleVelocities[i] = _propOfPic * _velocity(pos) + (1 - _propOfPic) * (_particleVelocities[i] + _deltaVelocity(pos));
	});
}

template <int Dim>
void ImplicitParticleLiquid<Dim>::maintainGridBasedData(StaggeredGridBasedScalarData<Dim> &weightSum)
{
	_deltaVelocity = _velocity;
	ParticleInCellLiquid<Dim>::maintainGridBasedData(weightSum);
}

template class ImplicitParticleLiquid<2>;
template class ImplicitParticleLiquid<3>;

}
