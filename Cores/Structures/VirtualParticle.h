#pragma once

#include "Structures/GridBasedScalarField.h"
#include "Structures/Particles.h"
#include "Structures/ParticlesBasedScalarField.h"
#include "Structures/ParticlesBasedVectorField.h"
#include "Structures/ParticlesNearbySearcher.h"
#include "Structures/SmoothedParticles.h"
#include "Structures/StaggeredGrid.h"

namespace PhysX {

    template<int Dim> class VirtualParticle : public SmoothedParticles<Dim> {
        DECLARE_DIM_TYPES(Dim)

    public:
        using SmoothedParticles<Dim>::positions;
        using SmoothedParticles<Dim>::densities;
        using SmoothedParticles<Dim>::volumes;

    protected:
        using SmoothedParticles<Dim>::_radius;
        using SmoothedParticles<Dim>::_kernelRadius;
        using SmoothedParticles<Dim>::_squaredKernelRadius;
        using SmoothedParticles<Dim>::_invKernelRadius;
        using SmoothedParticles<Dim>::_invSquaredKernelRadius;

        using SmoothedParticles<Dim>::_kernelNormCoeff0;
        using SmoothedParticles<Dim>::_kernelNormCoeff1;
        using SmoothedParticles<Dim>::_kernelNormCoeff2;

        using SmoothedParticles<Dim>::_nearbySearcher;

    private:
        using Particles<Dim>::clear;
        using Particles<Dim>::add;

    public:
        using SmoothedParticles<Dim>::resetNearbySearcher;
        using SmoothedParticles<Dim>::kernel;
        using SmoothedParticles<Dim>::gradientKernel;
        using SmoothedParticles<Dim>::parallelForEach;
        using SmoothedParticles<Dim>::forEachNearby;
        using SmoothedParticles<Dim>::forEach;
        using SmoothedParticles<Dim>::firstDerivativeKernel;

    private:
        const StaggeredGrid<Dim> _grid;
        GridBasedScalarData<Dim> _weight;

        real _alpha_0;
        real kappa;

        ParticlesBasedVectorField<Dim> _velocities;

        ParticlesBasedScalarField<Dim> _pressures;
        ParticlesBasedScalarField<Dim> _divergence;

        std::vector<Tripletr> _coefficients;
        SparseMatrixr         _matLaplacian;

    public:
        VirtualParticle(
            const StaggeredGrid<Dim> & grid,
            const real                 radius,
            const size_t               cnt = 0,
            const VectorDr &           pos = VectorDr::Zero());

        VirtualParticle & operator=(const VirtualParticle & rhs) = delete;
        virtual ~VirtualParticle()                               = default;

        virtual void computeDensities(const SmoothedParticles<Dim> & realParticles);
        virtual void computeVolumes(const SmoothedParticles<Dim> & realParticles);
        virtual void computeInfo(const SmoothedParticles<Dim> & realParticles);

        void setAlpha(const SmoothedParticles<Dim> & realParticles, const real target_rho);
        void setKappa(const real k) { kappa = k; }

        void generateParticles(const SmoothedParticles<Dim> & realParticles);

        void project(
            ParticlesBasedVectorField<Dim> &       velocity,
            const ParticlesBasedVectorField<Dim> & boundary_velocity,
            const SmoothedParticles<Dim> &         real_particle,
            const SmoothedParticles<Dim> &         boundary_particle,
            const real                             target_rho);

    protected:
        void buildLinearSystem(
            const ParticlesBasedVectorField<Dim> & velocity,
            const ParticlesBasedVectorField<Dim> & boundary_velocity,
            const SmoothedParticles<Dim> &         real_particle,
            const SmoothedParticles<Dim> &         boundary_particle,
            const real                             target_rho);
        void solveLinearSystem();
        void applyPressureGradient(
            ParticlesBasedVectorField<Dim> & velocity,
            const SmoothedParticles<Dim> &   real_particle,
            const real                       target_rho);

        void calculateVelocity(
            const ParticlesBasedVectorField<Dim> & velocity,
            const ParticlesBasedVectorField<Dim> & boundary_velocity,
            const SmoothedParticles<Dim> &         real_particle,
            const SmoothedParticles<Dim> &         boundary_particle);
        void calculateCoef(const real target_rho);
        void calculateDiv(
            const ParticlesBasedVectorField<Dim> & velocity,
            const ParticlesBasedVectorField<Dim> & boundary_velocity,
            const SmoothedParticles<Dim> &         real_particle,
            const SmoothedParticles<Dim> &         boundary_particle,
            const real                             target_rho);
    };

} // namespace PhysX
