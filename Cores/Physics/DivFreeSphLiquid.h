#pragma once

#include "Physics/SmthParticleHydrodLiquid.h"

namespace PhysX {

    template<int Dim> class DivFreeSphLiquid : public SmthParticleHydrodLiquid<Dim> {
        DECLARE_DIM_TYPES(Dim)

    protected:
        using SmthParticleHydrodLiquid<Dim>::_colliders;
        using SmthParticleHydrodLiquid<Dim>::_particles;
        using SmthParticleHydrodLiquid<Dim>::_velocities;
        using SmthParticleHydrodLiquid<Dim>::_pressures;
        using SmthParticleHydrodLiquid<Dim>::_enableGravity;
        using SmthParticleHydrodLiquid<Dim>::_viscosityCoeff;
        using SmthParticleHydrodLiquid<Dim>::_targetDensity;

    public:
        DivFreeSphLiquid(const real particleRadius): SmthParticleHydrodLiquid<Dim>(particleRadius) {}

        DivFreeSphLiquid(const DivFreeSphLiquid & rhs)             = delete;
        DivFreeSphLiquid & operator=(const DivFreeSphLiquid & rhs) = delete;
        virtual ~DivFreeSphLiquid()                                = default;

        virtual void initialize() override;

    protected:
        virtual void writeDescription(YAML::Node & root) const override;
        virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;

        void moveParticles(const real dt) override;

        void correctDensity();
        void correctVelocity();

        void updateVelocity(const real dt);
        void applyExternalForces(const real dt) override;
        void applyViscosityForce(const real dt) override;
        void applySurfaceTension(const real dt);
    };

} // namespace PhysX
