#include "SmoothedParticles.h"

#include "Structures/StaggeredGrid.h"

#include <numbers>

#include <cmath>

namespace PhysX {

    template<int Dim>
    SmoothedParticles<Dim>::SmoothedParticles(const real radius, const size_t cnt, const VectorDr & pos, const real mass):
        Particles<Dim>(cnt, pos, mass),
        _radius(radius),
        _kernelRadius(radius * 4),
        _squaredKernelRadius(_kernelRadius * _kernelRadius),
        _invKernelRadius(1 / _kernelRadius),
        _invSquaredKernelRadius(1 / _squaredKernelRadius),
        _kernelNormCoeff0(real(std::numbers::inv_pi) * _invSquaredKernelRadius * (Dim == 2 ? real(40) / 7 : 8 * _invKernelRadius)),
        _kernelNormCoeff1(6 * _invKernelRadius * _kernelNormCoeff0),
        _kernelNormCoeff2(2 * _invKernelRadius * _kernelNormCoeff1),
        _nearbySearcher(std::make_unique<HashGridSearcher<Dim>>(_kernelRadius)) {}

    template<int Dim>
    void SmoothedParticles<Dim>::computeDensities() {
        densities._data.resize(positions.size());
        parallelForEach([&](const int idx) {
            const VectorDr pos = positions[idx];
            densities[idx]     = 0;
            forEachNearby(pos, [&](const int j, const VectorDr & nearbyPos) {
                densities[idx] += kernel(nearbyPos - pos);
            });
            densities[idx] *= _mass;
        });
    }

    template<int Dim>
    void SmoothedParticles<Dim>::computeVolumes() {
        volumes._data.resize(positions.size());
        parallelForEach([&](const int idx) {
            volumes[idx] = _mass / densities[idx];
        });
    }

    template<int Dim>
    void SmoothedParticles<Dim>::computeVirtualVolumes(const SmoothedParticles & realParticles) {
        volumes._data.resize(positions.size());
        parallelForEach([&](const int I) {
            const VectorDr p_I = positions[I];
            volumes[I]         = 0;
            realParticles.forEachNearby(p_I, [&](const int j, const VectorDr & p_j) {
                volumes[I] += realParticles.volumes[j] * realParticles.volumes[j] * realParticles.kernel(p_I - p_j);
            });
        });
    }

    template<int Dim>
    void SmoothedParticles<Dim>::computeVirtualDensities(const SmoothedParticles & realParticles) {
        densities._data.resize(positions.size());
        parallelForEach([&](const int I) {
            const VectorDr p_I = positions[I];
            densities[I]       = 0;
            realParticles.forEachNearby(p_I, [&](const int j, const VectorDr & p_j) {
                densities[I] += kernel(p_I - p_j);
            });
            densities[I] *= realParticles._mass;
        });
    }

    template<int Dim>
    real SmoothedParticles<Dim>::getPackedKernelSum() const {
        if constexpr (Dim == 2)
            return kernel(0) + kernel(_radius * 2) * 6 + kernel(_radius * 2 * std::numbers::sqrt3) * 6;
        else
            return kernel(0) + kernel(_radius * 2) * 12 + kernel(_radius * 2 * std::numbers::sqrt2) * 6 + kernel(_radius * 2 * std::numbers::sqrt3) * 24;
    }

    template<int Dim>
    void SmoothedParticles<Dim>::generateBoxPacked(const VectorDr & center, const VectorDr & halfLengths) {
        if constexpr (Dim == 2) {
            const real     spacing = _radius * 2;
            const int      xScale  = int(std::round(halfLengths[0] * 2 / spacing));
            const int      yScale  = int(std::round((halfLengths[1] * 2 / spacing - 1) * 2 * real(std::numbers::inv_sqrt3) + 1));
            const VectorDr origin  = center - VectorDr(xScale - 1, (yScale - real(1)) * real(std::numbers::sqrt3) / 2) * spacing / 2;
            for (int j = 0; j < yScale; j++) {
                const VectorDr yOrigin = origin + VectorDr(j & 1 ? real(.5) : real(0), j * real(std::numbers::sqrt3) / 2) * spacing;
                for (int i = 0; i < xScale - (j & 1); i++)
                    positions._data.push_back(yOrigin + VectorDr::Unit(0) * i * spacing);
            }
        } else {
            const real         spacing    = 2 * real(std::numbers::sqrt2) * _radius;
            const VectorDi     resolution = ((halfLengths - VectorDr::Ones() * _radius) * 2 / spacing).array().round().template cast<int>().matrix();
            StaggeredGrid<Dim> grid(0, spacing, resolution, center);
            grid.forEachFace([&](const int axis, const VectorDi & face) {
                positions._data.push_back(grid.faceCenter(axis, face));
            });
            grid.forEachNode([&](const VectorDi & node) {
                positions._data.push_back(grid.nodeCenter(node));
            });
        }
    }

    template class SmoothedParticles<2>;
    template class SmoothedParticles<3>;

} // namespace PhysX
