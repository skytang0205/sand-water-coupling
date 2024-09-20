#pragma once

#include "Physics/SmthParticleHydrodLiquid.h"
#include "Structures/VirtualParticle.h"

namespace PhysX {

    template<int Dim> class DualParticleSphLiquid : public SmthParticleHydrodLiquid<Dim> {
        DECLARE_DIM_TYPES(Dim)

    protected:
        using SmthParticleHydrodLiquid<Dim>::_particles;
        using SmthParticleHydrodLiquid<Dim>::_velocities;
        using SmthParticleHydrodLiquid<Dim>::_colliders;

        VirtualParticle<Dim>           _virtual_particles;
        BoundaryParticles<Dim>         _boundary_particles;
        ParticlesBasedVectorField<Dim> _boundary_velocity;

        double _alpha_0;

        using SmthParticleHydrodLiquid<Dim>::_enableGravity;
        using SmthParticleHydrodLiquid<Dim>::_viscosityCoeff;
        using SmthParticleHydrodLiquid<Dim>::_targetDensity;

    public:
        DualParticleSphLiquid(const StaggeredGrid<Dim> & grid, const real particleRadius):
            SmthParticleHydrodLiquid<Dim>(particleRadius), _virtual_particles(grid, grid.spacing()),
            _boundary_particles(particleRadius), _alpha_0(-1e6) {}

        DualParticleSphLiquid(const DualParticleSphLiquid & rhs)             = delete;
        DualParticleSphLiquid & operator=(const DualParticleSphLiquid & rhs) = delete;
        virtual ~DualParticleSphLiquid()                                     = default;

        using SmthParticleHydrodLiquid<Dim>::getTimeStep;
        using SmthParticleHydrodLiquid<Dim>::dimension;

        virtual void writeDescription(YAML::Node & root) const override;
        virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;

        using SmthParticleHydrodLiquid<Dim>::saveFrame;
        using SmthParticleHydrodLiquid<Dim>::loadFrame;

        virtual void initialize() override;

        virtual void advance(const real dt) override;

    protected:
        using SmthParticleHydrodLiquid<Dim>::moveParticles;
        using SmthParticleHydrodLiquid<Dim>::applyExternalForces;
        using SmthParticleHydrodLiquid<Dim>::applyViscosityForce;

        virtual void applyPressureForce(const real dt);
    };

} // namespace PhysX
