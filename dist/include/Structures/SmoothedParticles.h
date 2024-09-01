#pragma once

#include "Geometries/Surface.h"
#include "Structures/Particles.h"
#include "Structures/ParticlesNearbySearcher.h"

namespace PhysX {

    template<int Dim> class SmoothedParticles : public Particles<Dim> {
        DECLARE_DIM_TYPES(Dim)

    public:
        using Particles<Dim>::positions;

        ParticlesScalarAttribute<Dim> densities;
        ParticlesScalarAttribute<Dim> volumes;

    private:
        using Particles<Dim>::_mass;
        using Particles<Dim>::_invMass;

    protected:
        const real _radius;
        const real _kernelRadius;
        const real _squaredKernelRadius;
        const real _invKernelRadius;
        const real _invSquaredKernelRadius;

        const real _kernelNormCoeff0;
        const real _kernelNormCoeff1;
        const real _kernelNormCoeff2;

        std::unique_ptr<ParticlesNearbySearcher<Dim>> _nearbySearcher;

    public:
        using Particles<Dim>::setMass;
        using Particles<Dim>::mass;
        using Particles<Dim>::forEach;
        using Particles<Dim>::parallelForEach;

        SmoothedParticles(
            const real radius, const size_t cnt = 0, const VectorDr & pos = VectorDr::Zero(), const real mass = 1);

        SmoothedParticles & operator=(const SmoothedParticles & rhs) = delete;
        virtual ~SmoothedParticles()                                 = default;

        real radius() const { return _radius; }
        real kernelRadius() const { return _kernelRadius; }

        virtual void computeDensities();
        virtual void computeVolumes();
        virtual void computeInfo();

        // Interpolation helper functions, using the cubic spline kernel

        real kernel(const real distance) const {
            const real x = distance * _invKernelRadius;
            if (x < real(.5)) return _kernelNormCoeff0 * (6 * x * x * (x - 1) + 1);
            else if (x < 1) return _kernelNormCoeff0 * 2 * (1 - x) * (1 - x) * (1 - x);
            else return 0;
        }

        real firstDerivativeKernel(const real distance) const {
            const real x = distance * _invKernelRadius;
            if (x < real(.5)) return _kernelNormCoeff1 * x * (3 * x - 2);
            else if (x < 1) return _kernelNormCoeff1 * (1 - x) * (x - 1);
            else return 0;
        }

        real secondDerivativeKernel(const real distance) const {
            const real x = distance * _invKernelRadius;
            if (x < real(.5)) return _kernelNormCoeff2 * (3 * x - 1);
            else if (x < 1) return _kernelNormCoeff2 * (1 - x);
            else return 0;
        }

        real laplacianKernel(const real distance) const {
            const real x = distance * _invKernelRadius;
            if (x < real(.5)) return _kernelNormCoeff2 * ((3 * x - 1) + (Dim - real(1)) / 2 * (3 * x - 2));
            else if (x < 1) return _kernelNormCoeff2 * ((1 - x) + (Dim - real(1)) / 2 / x * (1 - x) * (x - 1));
            else return 0;
        }

        real improvedLaplacianKernel(const real distance) const {
            const real x = distance * _invKernelRadius;
            if (x < real(.5)) return _kernelNormCoeff2 * (2 - 3 * x);
            else if (x < 1) return _kernelNormCoeff2 * (1 - x) * (1 - x) / x;
            else return 0;
        }

        real     kernel(const VectorDr & deltaPos) const { return kernel(deltaPos.norm()); }
        VectorDr gradientKernel(const VectorDr & deltaPos) const {
            return -firstDerivativeKernel(deltaPos.norm()) * deltaPos.normalized();
        }
        real laplacianKernel(const VectorDr & deltaPos) const { return laplacianKernel(deltaPos.norm()); }
        real improvedLaplacianKernel(const VectorDr & deltaPos) const {
            return improvedLaplacianKernel(deltaPos.norm());
        }

        real     getNeighborWeight(const VectorDr & pos) const;
        real     getBiasWeight(const VectorDr & pos) const;
        VectorDr getBiasGradient(const VectorDr & pos) const;

        real getPackedKernelSum() const;

        void resetNearbySearcher() { _nearbySearcher->reset(positions); }
        void forEachNearby(const VectorDr & pos, const std::function<void(const int, const VectorDr &)> & func) const {
            _nearbySearcher->forEach(positions, pos, func);
        }

        void generateBoxPacked(const VectorDr & center, const VectorDr & halfLengths);
    };

    template<int Dim> class BoundaryParticles : public SmoothedParticles<Dim> {
        DECLARE_DIM_TYPES(Dim)

    public:
        using SmoothedParticles<Dim>::positions;
        ParticlesVectorAttribute<Dim> norms;

    public:
        BoundaryParticles(
            const real radius, const size_t cnt = 0, const VectorDr & pos = VectorDr::Zero(), const real mass = 1):
            SmoothedParticles<Dim>(radius, cnt, pos, mass) {}

        BoundaryParticles & operator=(const BoundaryParticles & rhs) = delete;
        virtual ~BoundaryParticles()                                 = default;

        void addSurface(const Surface<Dim> & s, const VectorDr & center, const VectorDr & halfLengths);
    };

} // namespace PhysX
