#pragma once

#include "Geometries/Collider.h"
#include "Physics/Simulation.h"
#include "Structures/ParticlesBasedScalarField.h"
#include "Structures/ParticlesBasedVectorField.h"
#include "Structures/StaggeredGrid.h"
#include "Structures/VirtualParticle.h"
#include "Utilities/Shapes.h"

namespace PhysX {

    template<int Dim> class SmthParticleHydrodLiquid : public Simulation {
        DECLARE_DIM_TYPES(Dim)

    public:
        friend class SmthPartHydrodLiquidBuilder;

    protected:
        SmoothedParticles<Dim>         _particles;
        ParticlesBasedVectorField<Dim> _velocities;

        double _alpha_0;

        VirtualParticle<Dim> _virtual_particles;

        BoundaryParticles<Dim>         _boundary_particles;
        ParticlesBasedVectorField<Dim> _boundary_velocity;

        std::vector<std::unique_ptr<Collider<Dim>>> _colliders;

        bool _enableGravity = false;

        real _viscosityCoeff = 0;
        // real _viscosityCoeff = real(.01);
        real _targetDensity = real(1e3);
        real _eosMultiplier = 25;

    public:
        SmthParticleHydrodLiquid(const StaggeredGrid<Dim> & grid, const real particleRadius):
            _particles(particleRadius), _virtual_particles(grid, grid.spacing()), _boundary_particles(particleRadius),
            _alpha_0(-1e6) {}

        SmthParticleHydrodLiquid(const SmthParticleHydrodLiquid & rhs)             = delete;
        SmthParticleHydrodLiquid & operator=(const SmthParticleHydrodLiquid & rhs) = delete;
        virtual ~SmthParticleHydrodLiquid()                                        = default;

        virtual real getTimeStep(const uint frameRate, const real stepRate) const override {
            return stepRate * _particles.radius() * 2 / _velocities.normMax();
        }

        virtual int  dimension() const override { return Dim; }
        virtual void writeDescription(YAML::Node & root) const override;
        virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;
        virtual void saveFrame(const std::string & frameDir) const override;
        virtual void loadFrame(const std::string & frameDir) override;

        virtual void initialize() override;
        virtual void advance(const real dt) override;

        void generateSurface(const Surface<Dim> & surface);
        void addShape(const Shapes<Dim> & shape);

    protected:
        virtual void reinitializeParticlesBasedData();
        virtual void moveParticles(const real dt);
        virtual void applyExternalForces(const real dt);
        virtual void applyViscosityForce(const real dt);
        virtual void applyPressureForce(const real dt);
    };

} // namespace PhysX
