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
	{ // Save markerVelocities.
		std::ofstream fout(frameDir + "/markerVelocities.sav", std::ios::binary);
		_markerVelocities.save(fout);
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
	{ // Load markerVelocities.
		std::ifstream fin(frameDir + "/markerVelocities.sav", std::ios::binary);
		_markerVelocities.load(fin);
	}
	{ // Load delteVelocity.
		std::ifstream fin(frameDir + "/delteVelocity.sav", std::ios::binary);
		_deltaVelocity.load(fin);
	}
}

template <int Dim>
void ImplicitParticleLiquid<Dim>::transferFromGridsToParticles()
{
	_deltaVelocity.parallelForEach([&](const int axis, const VectorDi &face) {
		_deltaVelocity[axis][face] = _velocity[axis][face] - _deltaVelocity[axis][face];
	});
	_markerVelocities.parallelForEach([&](const int i) {
		const VectorDr pos = _markerPositions[i];
		_markerVelocities[i] = _propOfPic * _velocity(pos) + (1 - _propOfPic) * (_markerVelocities[i] + _deltaVelocity(pos));
	});
}

template <int Dim>
void ImplicitParticleLiquid<Dim>::maintainGridsBasedData(StaggeredGridBasedData<Dim> &weightSum)
{
	_deltaVelocity = _velocity;
	ParticleInCellLiquid<Dim>::maintainGridsBasedData(weightSum);
}

template class ImplicitParticleLiquid<2>;
template class ImplicitParticleLiquid<3>;

}
