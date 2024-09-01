#pragma once

#include "Geometries/Collider.h"
#include "Physics/Simulation.h"
#include "Structures/ParticlesBasedScalarField.h"
#include "Structures/ParticlesBasedVectorField.h"
#include "Structures/StaggeredGrid.h"
#include "Structures/DEMParticle.h"

namespace PhysX {

    template<int Dim> class DEMParticleSand : public Simulation {
        DECLARE_DIM_TYPES(Dim)

    public:
        friend class DEMParticleSandBuilder;

    protected:
        DEMParticle<Dim>         _particles;

        BoundaryParticles<Dim>         _boundary_particles;
        ParticlesBasedVectorField<Dim> _boundary_velocity;

        std::vector<std::unique_ptr<Collider<Dim>>> _colliders;

        bool _enableGravity = false;

        real _targetDensity = real(1e3);
        real _eosMultiplier = 25;

    public:
        DEMParticleSand(const real particleRadius):
            _particles(particleRadius), _boundary_particles(particleRadius),
            _alpha_0(-1e6) {}

        DEMParticleSandd(const DEMParticleSand & rhs)             = delete;
        DEMParticleSand & operator=(const DEMParticleSand & rhs) = delete;
        virtual ~DEMParticleSand()                                        = default;

        virtual real getTimeStep(const uint frameRate, const real stepRate) const override {
            return stepRate * _particles.radius() * 2 / _particles.velosities.normMax();
        }

        virtual int  dimension() const override { return Dim; }
        virtual void writeDescription(YAML::Node & root) const override;
        virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;
        virtual void saveFrame(const std::string & frameDir) const override;
        virtual void loadFrame(const std::string & frameDir) override;

        virtual void initialize() override;
        virtual void advance(const real dt) override;

        void generateSurface(const Surface<Dim> & surface);

    protected:
        virtual void reinitializeParticlesBasedData();
        virtual void moveParticles(const real dt);
        virtual void applyExternalForces(const real dt);
        virtual void applyViscosityForce(const real dt);
        virtual void applyPressureForce(const real dt);
    };

} // namespace PhysX
