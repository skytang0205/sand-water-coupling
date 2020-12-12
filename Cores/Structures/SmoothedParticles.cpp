#include "SmoothedParticles.h"

#include "Structures/StaggeredGrid.h"

#include <numbers>

#include <cmath>

namespace PhysX {

template <int Dim>
SmoothedParticles<Dim>::SmoothedParticles(const real radius, const size_t cnt, const VectorDr &pos, const real mass) :
	Particles<Dim>(cnt, pos, mass),
	_radius(radius),
	_kernelRadius(radius * 4),
	_squaredKernelRadius(_kernelRadius *_kernelRadius),
	_invKernelRadius(1 / _kernelRadius),
	_invSquaredKernelRadius(1 / _squaredKernelRadius),
	_stdKernelNormCoeff0((Dim == 2 ? real(4) : 315 * _invKernelRadius / 64) * real(std::numbers::inv_pi) * _invSquaredKernelRadius),
	_stdKernelNormCoeff1(-6 * _invSquaredKernelRadius * _stdKernelNormCoeff0),
	_stdKernelNormCoeff2(_stdKernelNormCoeff1),
	_spikyKernelNormCoeff0((Dim == 2 ? real(10) : 15 * _invKernelRadius) * real(std::numbers::inv_pi) *_invSquaredKernelRadius),
	_spikyKernelNormCoeff1(-3 * _invKernelRadius * _spikyKernelNormCoeff0),
	_spikyKernelNormCoeff2(-2 * _invKernelRadius * _spikyKernelNormCoeff1),
	_nearbySearcher(std::make_unique<HashGridSearcher<Dim>>(_kernelRadius))
{ }

template <int Dim>
void SmoothedParticles<Dim>::computeDensities()
{
	densities._data.resize(positions.size());
	parallelForEach([&](const int idx) {
		const VectorDr pos = positions[idx];
		densities[idx] = 0;
		forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
			densities[idx] += stdKernel(nearbyPos - pos);
		});
		densities[idx] *= _mass;
	});
}

template <int Dim>
real SmoothedParticles<Dim>::getPackedStdKernelSum() const
{
	if constexpr (Dim == 2)
		return stdKernel(0) + stdKernel(_radius * 2) * 6 + stdKernel(_radius * 2 * std::numbers::sqrt3) * 6;
	else
		return 0;
}

template <int Dim>
void SmoothedParticles<Dim>::generateBoxPacked(const VectorDr &center, const VectorDr &halfLengths)
{
	const real spacing = _radius * 2;
	if constexpr (Dim == 2) {
		const int xScale = int(std::round(halfLengths[0] * 2 / spacing));
		const int yScale = int(std::round((halfLengths[1] * 2 / spacing - 1) * 2 * real(std::numbers::inv_sqrt3) + 1));
		const VectorDr origin = -VectorDr(xScale - 1, (yScale - real(1)) * real(std::numbers::sqrt3) / 2) * spacing / 2;
		for (int j = 0; j < yScale; j++) {
			const VectorDr yOrigin = origin + VectorDr(j & 1 ? real(.5) : real(0), j * real(std::numbers::sqrt3) / 2) * spacing;
			for (int i = 0; i < xScale - (j & 1); i++)
				positions._data.push_back(yOrigin + VectorDr::Unit(0) * i * spacing);
		}
	}
	else {
	}
}

template class SmoothedParticles<2>;
template class SmoothedParticles<3>;

}
