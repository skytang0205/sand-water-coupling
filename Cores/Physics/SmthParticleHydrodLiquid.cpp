#include "SmthParticleHydrodLiquid.h"

#include "Utilities/Constants.h"

namespace PhysX {

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::writeDescription(YAML::Node & root) const {
        { // Description of particles.
            YAML::Node node;
            node["name"]                       = "particles";
            node["data_mode"]                  = "dynamic";
            node["primitive_type"]             = "point_list";
            node["material"]["diffuse_albedo"] = (Vector4f(52, 108, 156, 255) / 255).eval(); // Haijun Blue
            node["indexed"]                    = false;
            node["color_map"]["enabled"]       = true;
            root["objects"].push_back(node);
        }
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::writeFrame(const std::string & frameDir, const bool staticDraw) const {
        { // Write particles.
            std::ofstream fout(frameDir + "/particles.mesh", std::ios::binary);
            IO::writeValue(fout, uint(_particles.size()));
            _particles.forEach(
                [&](const int i) { IO::writeValue(fout, _particles.positions[i].template cast<float>().eval()); });
            if constexpr (Dim == 3) {
                _particles.forEach([&](const int i) { IO::writeValue(fout, VectorDf::Unit(2).eval()); });
            }
            _particles.forEach([&](const int i) { IO::writeValue(fout, float(_velocities[i].norm())); });
        }
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::saveFrame(const std::string & frameDir) const {
        { // Save particles.
            std::ofstream fout(frameDir + "/particles.sav", std::ios::binary);
            _particles.positions.save(fout);
        }
        { // save velocities.
            std::ofstream fout(frameDir + "/velocities.sav", std::ios::binary);
            _velocities.save(fout);
        }
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::loadFrame(const std::string & frameDir) {
        { // Load particles.
            std::ifstream fin(frameDir + "/particles.sav", std::ios::binary);
            _particles.positions.load(fin);
        }
        reinitializeParticlesBasedData();
        { // Load velocities.
            std::ifstream fin(frameDir + "/velocities.sav", std::ios::binary);
            _velocities.load(fin);
        }
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::initialize() { reinitializeParticlesBasedData(); }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::advance(const real dt) {
        moveParticles(dt);

        applyExternalForces(dt);
        applyViscosityForce(dt);

        applyPressureForce(dt);
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::reinitializeParticlesBasedData() {
        _pressures.resize(&_particles);
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::moveParticles(const real dt) {
        _particles.parallelForEach([&](const int i) { _particles.positions[i] += _velocities[i] * dt; });

        // Resolve collisions.
        for (const auto & collider : _colliders) {
            collider->collide(_particles.positions, _velocities, _particles.radius());
        }

        _particles.resetNearbySearcher();
        _particles.computeDensities();
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::applyExternalForces(const real dt) {
        if (_enableGravity) {
            _particles.parallelForEach([&](const int i) { _velocities[i][1] -= kGravity * dt; });
        }
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::applyViscosityForce(const real dt) {
        if (_viscosityCoeff) {
            auto newVelocities = _velocities;
            _particles.parallelForEach([&](const int i) {
                newVelocities[i] += _viscosityCoeff * _velocities.divFreeLaplacianAtDataPoint(i) * dt;
            });
            _velocities = newVelocities;
        }
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::applyPressureForce(const real dt) {
        // Compute pressures by the equation of state.
        _particles.parallelForEach([&](const int i) {
            _pressures[i] = _eosMultiplier * (_particles.densities[i] - _targetDensity);
            if (_pressures[i] < 0) _pressures[i] = 0;
        });

        // Apply pressure gradient.
        _particles.parallelForEach([&](const int i) {
            _velocities[i] -= _pressures.symmetricGradientAtDataPoint(i) / _particles.densities[i] * dt;
        });
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::generateSurface(const Surface<Dim> & surface) {
        //_virtual_weight.forEach([&](const VectorDi & nodeCell) {
        //    VectorDr pos = _virtual_weight.position(nodeCell);
        //    if (surface.isInside(pos)) _boundary_particles.positions._data.push_back(pos);
        //});

        //_boundary_particles._mass = _particles._mass;
        //_boundary_velocity.resize(&_boundary_particles);
        //_boundary_particles.resetNearbySearcher();
        //_boundary_particles.computeDensities();
        //_boundary_particles.computeVolumes();
    }

    template<int Dim> void SmthParticleHydrodLiquid<Dim>::addShape(const Shapes<Dim> & shape) {
        auto origin_size = _particles.positions.size();
        _particles.resize(origin_size + shape.positions.size());
        _velocities.resize(&_particles);
        for (int i = 0; i < shape.positions.size(); i++) {
            _particles.positions[origin_size + i] = shape.positions[i];
            _velocities[origin_size + i]          = shape.velocities[i];
        }
    }

    template class SmthParticleHydrodLiquid<2>;
    template class SmthParticleHydrodLiquid<3>;

} // namespace PhysX
