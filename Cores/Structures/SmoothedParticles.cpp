#include "SmoothedParticles.h"

#include "Structures/StaggeredGrid.h"

#include "Solvers/IterativeSolver.h"

#include <numbers>

#include <cmath>

namespace PhysX {

    template<int Dim>
    SmoothedParticles<Dim>::SmoothedParticles(
        const real radius, const size_t cnt, const VectorDr & pos, const real mass, const real kernel_rate):
        Particles<Dim>(cnt, pos, mass),
        _radius(radius), _kernelRadius(radius * kernel_rate), _squaredKernelRadius(_kernelRadius * _kernelRadius),
        _invKernelRadius(1 / _kernelRadius), _invSquaredKernelRadius(1 / _squaredKernelRadius),
        _kernelNormCoeff0(
            real(std::numbers::inv_pi) * _invSquaredKernelRadius * (Dim == 2 ? real(40) / 7 : 8 * _invKernelRadius)),
        _kernelNormCoeff1(6 * _invKernelRadius * _kernelNormCoeff0),
        _kernelNormCoeff2(2 * _invKernelRadius * _kernelNormCoeff1),
        _nearbySearcher(std::make_unique<HashGridSearcher<Dim>>(_kernelRadius)) {}

    template<int Dim> void SmoothedParticles<Dim>::computeDensities() {
        densities._data.resize(positions.size());
        parallelForEach([&](const int idx) {
            const VectorDr pos = positions[idx];
            densities[idx]     = _mass * getNeighborWeight(pos);
        });
    }

    template<int Dim> void SmoothedParticles<Dim>::computeVolumes() {
        volumes._data.resize(positions.size());
        parallelForEach([&](const int idx) { volumes[idx] = _mass / densities[idx]; });
    }

    template<int Dim> void SmoothedParticles<Dim>::computeInfo() {
        computeDensities();
        computeVolumes();
    }

    template<int Dim> void SmoothedParticles<Dim>::resize(const size_t cnt, const VectorDr & pos) {
        positions._data.resize(cnt, pos);
        densities._data.resize(cnt, 0);
    }

    template<int Dim> real SmoothedParticles<Dim>::getNeighborWeight(const VectorDr & pos) const {
        double sum = 0;

        forEachNearby(pos, [&](const int j, const VectorDr & nearbyPos) { sum += kernel(pos - nearbyPos); });

        return sum;
    }

    template<int Dim> real SmoothedParticles<Dim>::getBiasWeight(const VectorDr & pos) const {
        double sum = 0;

        forEachNearby(
            pos, [&](const int j, const VectorDr & nearbyPos) { sum += volumes[j] * kernel(pos - nearbyPos); });

        return sum;
    }

    template<int Dim> auto SmoothedParticles<Dim>::getBiasGradient(const VectorDr & pos) const -> VectorDr {
        VectorDr v = VectorDr::Zero();

        forEachNearby(
            pos, [&](const int j, const VectorDr & nearbyPos) { v += volumes[j] * gradientKernel(pos - nearbyPos); });

        return v;
    }

    template<int Dim> real SmoothedParticles<Dim>::getPackedKernelSum() const {
        if constexpr (Dim == 2)
            return kernel(0) + kernel(_radius * 2) * 6 + kernel(_radius * 2 * std::numbers::sqrt3) * 6;
        else
            return kernel(0) + kernel(_radius * 2) * 12 + kernel(_radius * 2 * std::numbers::sqrt2) * 6
                + kernel(_radius * 2 * std::numbers::sqrt3) * 24;
    }

    template<int Dim>
    void SmoothedParticles<Dim>::generateBoxPacked(const VectorDr & center, const VectorDr & halfLengths) {
        if constexpr (Dim == 2) {
            const real spacing = _radius * 2;
            const int  xScale  = int(std::round(halfLengths[0] * 2 / spacing));
            const int  yScale =
                int(std::round((halfLengths[1] * 2 / spacing - 1) * 2 * real(std::numbers::inv_sqrt3) + 1));
            const VectorDr origin =
                center - VectorDr(xScale - 1, (yScale - real(1)) * real(std::numbers::sqrt3) / 2) * spacing / 2;
            for (int j = 0; j < yScale; j++) {
                const VectorDr yOrigin =
                    origin + VectorDr(j & 1 ? real(.5) : real(0), j * real(std::numbers::sqrt3) / 2) * spacing;
                for (int i = 0; i < xScale - (j & 1); i++)
                    positions._data.push_back(yOrigin + VectorDr::Unit(0) * i * spacing);
            }
        } else {
            const real     spacing    = 2 * real(std::numbers::sqrt2) * _radius;
            const VectorDi resolution = ((halfLengths - VectorDr::Ones() * _radius) * 2 / spacing)
                                            .array()
                                            .round()
                                            .template cast<int>()
                                            .matrix();
            StaggeredGrid<Dim> grid(0, spacing, resolution, center);
            grid.forEachFace(
                [&](const int axis, const VectorDi & face) { positions._data.push_back(grid.faceCenter(axis, face)); });
            grid.forEachNode([&](const VectorDi & node) { positions._data.push_back(grid.nodeCenter(node)); });
        }
    }

    template<int Dim>
    void BoundaryParticles<Dim>::addSurface(
        const Surface<Dim> & s, const VectorDr & center, const VectorDr & halfLengths) {
        if constexpr (Dim == 2) {
            const real spacing = this->_radius * 2;
            const int  xScale  = int(std::round(halfLengths[0] * 2 / spacing));
            const int  yScale =
                int(std::round((halfLengths[1] * 2 / spacing - 1) * 2 * real(std::numbers::inv_sqrt3) + 1));
            const VectorDr origin =
                center - VectorDr(xScale - 1, (yScale - real(1)) * real(std::numbers::sqrt3) / 2) * spacing / 2;
            for (int j = 0; j < yScale; j++) {
                const VectorDr yOrigin =
                    origin + VectorDr(j & 1 ? real(.5) : real(0), j * real(std::numbers::sqrt3) / 2) * spacing;
                for (int i = 0; i < xScale - (j & 1); i++) {
                    VectorDr pos = yOrigin + VectorDr::Unit(0) * i * spacing;
                    if (s.isInside(pos)) {
                        positions._data.push_back(pos);
                        norms._data.push_back(s.closestNormal(pos));
                    }
                }
            }
        } else {
            const real     spacing    = 2 * real(std::numbers::sqrt2) * this->_radius;
            const VectorDi resolution = ((halfLengths - VectorDr::Ones() * this->_radius) * 2 / spacing)
                                            .array()
                                            .round()
                                            .template cast<int>()
                                            .matrix();
            StaggeredGrid<Dim> grid(0, spacing, resolution, center);
            grid.forEachFace([&](const int axis, const VectorDi & face) {
                VectorDr pos = grid.faceCenter(axis, face);
                if (s.isInside(pos)) {
                    positions._data.push_back(pos);
                    norms._data.push_back(s.closestNormal(pos));
                }
            });
            grid.forEachNode([&](const VectorDi & node) {
                VectorDr pos = grid.nodeCenter(node);
                if (s.isInside(pos)) {
                    positions._data.push_back(pos);
                    norms._data.push_back(s.closestNormal(pos));
                }
            });
        }
    }

    template<int Dim>
    WeakCompParticles<Dim>::WeakCompParticles(
        const real radius, const size_t cnt, const VectorDr & pos, const real mass, const real kernel_rate):
        SmoothedParticles<Dim>(radius, cnt, pos, mass),
        _kernelNormCoeff(
            real(std::numbers::inv_pi) * _invSquaredKernelRadius * (Dim == 2 ? 10. / 7. : _invKernelRadius)) {}

    template class SmoothedParticles<2>;
    template class SmoothedParticles<3>;
    template class BoundaryParticles<2>;
    template class BoundaryParticles<3>;
    template class WeakCompParticles<2>;
    template class WeakCompParticles<3>;

} // namespace PhysX
