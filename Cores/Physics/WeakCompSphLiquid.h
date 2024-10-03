#pragma once

#include "Physics/SmthParticleHydrodLiquid.h"

namespace PhysX {

    template<int Dim> class WeakCompSphLiquid : public SmthParticleHydrodLiquid<Dim> {
        DECLARE_DIM_TYPES(Dim)

    protected:
        using SmthParticleHydrodLiquid<Dim>::_colliders;
        using SmthParticleHydrodLiquid<Dim>::_velocities;
        using SmthParticleHydrodLiquid<Dim>::_pressures;
        using SmthParticleHydrodLiquid<Dim>::_enableGravity;
        using SmthParticleHydrodLiquid<Dim>::_viscosityCoeff;
        using SmthParticleHydrodLiquid<Dim>::_targetDensity;

        WeakCompParticles<Dim> _particles;

        real _adiabaticIndex = 7.;
        real _c_s;
        real _pressureCoeff;
        real _surfaceTension = 0.01;

    public:
        WeakCompSphLiquid(const real particleRadius):
            SmthParticleHydrodLiquid<Dim>(particleRadius), _particles(particleRadius), _c_s(5),
            _pressureCoeff(1e3 * _c_s * _c_s / _adiabaticIndex) {}

        WeakCompSphLiquid(const WeakCompSphLiquid & rhs)             = delete;
        WeakCompSphLiquid & operator=(const WeakCompSphLiquid & rhs) = delete;
        virtual ~WeakCompSphLiquid()                                 = default;

        virtual void initialize() override;

        virtual void advance(const real dt) override;

        virtual void addShape(const Shapes<Dim> & shape);

    protected:
        virtual void writeDescription(YAML::Node & root) const override;
        virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;

        void moveParticles(const real dt) override;

        void calculatePressure();

        void updateVelocity(const real dt);
        void applyExternalForces(const real dt) override;
        void applyViscosityForce(const real dt) override;
        void applySurfaceTension(const real dt);
        void applyPressure(const real dt);

        void updateDensity(const real dt);
    };

} // namespace PhysX
