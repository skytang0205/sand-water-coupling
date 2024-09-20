#include "DEMParticleSand.h"

#include "Utilities/Constants.h"

namespace PhysX {

    template<int Dim> void DEMParticleSand<Dim>::writeDescription(YAML::Node & root) const {
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
        { // Description of particles.
            YAML::Node node;
            node["name"]                       = "boundary_particles";
            node["data_mode"]                  = "dynamic";
            node["primitive_type"]             = "point_list";
            node["material"]["diffuse_albedo"] = (Vector4f(52, 108, 156, 255) / 255).eval(); // Haijun Blue
            node["indexed"]                    = false;
            node["color_map"]["enabled"]       = false;
            root["objects"].push_back(node);
        }
    }

    template<int Dim>
    void DEMParticleSand<Dim>::writeFrame(const std::string & frameDir, const bool staticDraw) const {
        { // Write particles.
            std::ofstream fout(frameDir + "/particles.mesh", std::ios::binary);
            IO::writeValue(fout, uint(_particles.size()));
            _particles.forEach(
                [&](const int i) { IO::writeValue(fout, _particles.positions[i].template cast<float>().eval()); });
            if constexpr (Dim == 3) {
                _particles.forEach([&](const int i) { IO::writeValue(fout, VectorDf::Unit(2).eval()); });
            }
            _particles.forEach([&](const int i) { IO::writeValue(fout, float(_particles.velocities[i].norm())); });
        }
        { // Write particles.
            std::ofstream fout(frameDir + "/boundary_particles.mesh", std::ios::binary);
            IO::writeValue(fout, uint(_boundary_particles.size()));
            _boundary_particles.forEach([&](const int i) {
                IO::writeValue(fout, _boundary_particles.positions[i].template cast<float>().eval());
            });
            if constexpr (Dim == 3) {
                _boundary_particles.forEach([&](const int i) { IO::writeValue(fout, VectorDf::Unit(2).eval()); });
            }
            //_particles.forEach([&](const int i) { IO::writeValue(fout, float(_velocities[i].norm())); });
        }
    }

    template<int Dim> void DEMParticleSand<Dim>::saveFrame(const std::string & frameDir) const {
        { // Save particles.
            std::ofstream fout(frameDir + "/particles.sav", std::ios::binary);
            _particles.positions.save(fout);
        }
        { // save velocities.
            std::ofstream fout(frameDir + "/velocities.sav", std::ios::binary);
            _particles.velocities.save(fout);
        }
    }

    template<int Dim> void DEMParticleSand<Dim>::loadFrame(const std::string & frameDir) {
        { // Load particles.
            std::ifstream fin(frameDir + "/particles.sav", std::ios::binary);
            _particles.positions.load(fin);
        }
        reinitializeParticlesBasedData();
        { // Load velocities.
            std::ifstream fin(frameDir + "/velocities.sav", std::ios::binary);
            _particles.velocities.load(fin);
        }
    }

    template<int Dim> void DEMParticleSand<Dim>::initialize() {
        _boundary_particles.setMass(_particles.mass());

        _particles.resetNearbySearcher();
        _particles.computeInfo();

        _boundary_particles.resetNearbySearcher();
        _boundary_particles.computeInfo();

        std::cout << "Boundary Particle size: " << _boundary_particles.size() << std::endl;
    }

    template<int Dim> void DEMParticleSand<Dim>::advance(const real dt) {
        moveParticles(dt);
        applyPressureForce(dt);
        applyExternalForces(dt);        
    }

    template<int Dim> void DEMParticleSand<Dim>::reinitializeParticlesBasedData() {}

    template<int Dim> void DEMParticleSand<Dim>::moveParticles(const real dt) {
        _particles.parallelForEach([&](const int i) { _particles.positions[i] += _particles.velocities[i] * dt; });

        // Resolve collisions.
        // for (const auto & collider : _colliders) {
        //    collider->collide(_particles.positions, _velocities, _particles.radius());
        //}

        _particles.resetNearbySearcher();
        //_particles.computeDensities();
        //_particles.computeVolumes();
    }

    template<int Dim> void DEMParticleSand<Dim>::applyExternalForces(const real dt) {
        if (_enableGravity) {
            _particles.parallelForEach([&](const int i) { _particles.velocities[i][1] -= kGravity * dt; });
        }
    }


    template<int Dim> void DEMParticleSand<Dim>::applyPressureForce(const real dt) {
        //// Compute pressures by the equation of state.
        //_particles.parallelForEach([&](const int i) {
        //    _pressures[i] = _eosMultiplier * (_particles.densities[i] -
        //    _targetDensity); if (_pressures[i] < 0) _pressures[i] = 0;
        //});

        //// Apply pressure gradient.
        //_particles.parallelForEach([&](const int i) {
        //    _velocities[i] -= _pressures.symmetricGradientAtDataPoint(i) /
        //    _particles.densities[i] * dt;
        //});
        _particles.parallelForEach([&](const int i) {
            VectorDr p_i   = _particles.positions[i];
            VectorDr force = VectorDr::Zero();
            force += _particles.getForceSum(i) * dt;
            //std::cout << getForceSum(i) << std::endl;

            //force = VectorDr::Zero(); 

            _boundary_particles.forEachNearby(p_i, [&](const int js, const VectorDr & p_js) {
                VectorDr dis_vel = _particles.velocities[i] - _boundary_velocity[js];
                real dis_vel_norm = _boundary_particles.norms[js].dot(dis_vel);
                VectorDr F_norm = - _boundary_particles.volumes[js]
                    * std::min(dis_vel_norm, 0.)
                    * _boundary_particles.norms[js] * _boundary_particles.kernel(p_i - p_js);
                VectorDr F_tang = - _particles._K_tang() * (dis_vel - dis_vel_norm * _boundary_particles.norms[js]);
                if(F_tang.norm() <= 0.000001)
                    force += F_norm + F_tang;
                else
                    force += F_norm + (F_norm.norm() * _particles._tan_fricangle() > F_tang.norm() ?  F_tang : (F_norm.norm() * _particles._tan_fricangle() / F_tang.norm()) * F_tang);
                });
            _particles.velocities[i] += force;
        });
    }

    template<int Dim> void DEMParticleSand<Dim>::generateSurface(const Surface<Dim> & surface) {
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

    template<int Dim> void DEMParticleSand<Dim>::addShape(const Shapes<Dim> & shape) {
        auto origin_size =  _particles.positions.size();
        _particles.resize(origin_size + shape.positions.size());
        _particles.velocities.resize(&_particles);
        for (int i = 0; i < shape.positions.size(); i++){
            _particles.positions[origin_size + i] = shape.positions[i];
            _particles.velocities[origin_size + i] = shape.velocities[i];
        }
    } 

    template class DEMParticleSand<2>;
    template class DEMParticleSand<3>;

} // namespace PhysX
