#pragma once

#include "Geometries/Collider.h"
#include "Physics/Simulation.h"
#include "Structures/GridBasedScalarField.h"
#include "Structures/ParticlesBasedScalarField.h"
#include "Structures/ParticlesBasedVectorField.h"
#include "Structures/StaggeredGrid.h"

namespace PhysX {

    template<int Dim>
    class SmthParticleHydrodLiquid : public Simulation {
        DECLARE_DIM_TYPES(Dim)

    public:
        friend class SmthPartHydrodLiquidBuilder;

    protected:
        SmoothedParticles<Dim>         _particles;
        ParticlesBasedVectorField<Dim> _velocities;

        const StaggeredGrid<Dim> _grid;
        GridBasedScalarData<Dim> _virtual_weight;

        double _alpha_0;

        SmoothedParticles<Dim>         _virtual_particles;
        ParticlesBasedScalarField<Dim> _pressures;
        ParticlesBasedScalarField<Dim> _divergence;
        ParticlesBasedVectorField<Dim> _virtual_velocity;

        std::vector<Tripletr> _coefficients;
        SparseMatrixr         _matLaplacian;

        std::vector<std::unique_ptr<Collider<Dim>>> _colliders;

        bool _enableGravity = false;

        real _viscosityCoeff = 0;
        // real _viscosityCoeff = real(.01);
        real _targetDensity = real(1e3);
        real _eosMultiplier = 25;

    public:
        SmthParticleHydrodLiquid(const StaggeredGrid<Dim> & grid, const real particleRadius):
            _grid(grid),
            _virtual_weight(_grid.nodeGrid(), 0),
            _particles(particleRadius),
            _virtual_particles(_grid.spacing()),
            _alpha_0(-1e6) {}

        SmthParticleHydrodLiquid(const SmthParticleHydrodLiquid & rhs)             = delete;
        SmthParticleHydrodLiquid & operator=(const SmthParticleHydrodLiquid & rhs) = delete;
        virtual ~SmthParticleHydrodLiquid()                                        = default;

        virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return stepRate * _particles.radius() * 2 / _velocities.normMax(); }

        virtual int  dimension() const override { return Dim; }
        virtual void writeDescription(YAML::Node & root) const override;
        virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;
        virtual void saveFrame(const std::string & frameDir) const override;
        virtual void loadFrame(const std::string & frameDir) override;

        virtual void initialize() override;
        virtual void advance(const real dt) override;

    protected:
        virtual void reinitializeParticlesBasedData();
        virtual void moveParticles(const real dt);
        virtual void applyExternalForces(const real dt);
        virtual void applyViscosityForce(const real dt);
        virtual void applyPressureForce(const real dt);

        void generateVirtualParticle();

        void calculateCoef();
        void calculateVirtualVelocity();
        void calculateDiv(real dt);
    };

} // namespace PhysX
