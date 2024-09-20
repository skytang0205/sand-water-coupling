#include "WeakCompSphLiquid.h"

#include "Utilities/Constants.h"

namespace PhysX {

    template<int Dim> void WeakCompSphLiquid<Dim>::writeDescription(YAML::Node & root) const {
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
    void WeakCompSphLiquid<Dim>::writeFrame(const std::string & frameDir, const bool staticDraw) const {
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

    template<int Dim> inline void WeakCompSphLiquid<Dim>::initialize() {
        _pressures.resize(&_particles);

        _particles.parallelForEach([&](const int i) { _particles.densities[i] = _targetDensity; });
        _velocities.setZero();
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::addShape(const Shapes<Dim> & shape) {
        auto origin_size = _particles.positions.size();
        _particles.resize(origin_size + shape.positions.size());
        _velocities.resize(&_particles);
        for (int i = 0; i < shape.positions.size(); i++) {
            _particles.positions[origin_size + i] = shape.positions[i];
            _velocities[origin_size + i]          = shape.velocities[i];
        }
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::advance(const real dt) {
        moveParticles(dt);

        calculatePressure();

        updateVelocity(dt);

        updateDensity(dt);
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::moveParticles(const real dt) {
        _particles.parallelForEach([&](const int i) { _particles.positions[i] += _velocities[i] * dt; });

        // Resolve collisions.
        for (const auto & collider : _colliders) {
            collider->collide(_particles.positions, _velocities, _particles.radius());
        }

        _particles.resetNearbySearcher();
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::calculatePressure() {
        _particles.parallelForEach([&](const int i) {
            _pressures[i] = _pressureCoeff * (std::pow(_particles.densities[i] / _targetDensity, _adiabaticIndex) - 1.);
        });
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::updateVelocity(const real dt) {
        applyPressure(dt);
        applyExternalForces(dt);
        applyViscosityForce(dt);
        applySurfaceTension(dt);
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::applyExternalForces(const real dt) {
        if (_enableGravity) {
            _particles.parallelForEach([&](const int i) { _velocities[i][1] -= kGravity * dt; });
        }
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::applyViscosityForce(const real dt) {
        if (_viscosityCoeff) {
            double h = _particles.kernelRadius();

            ParticlesBasedVectorData<Dim> _deltaVelocity;
            _deltaVelocity.resize(&_particles);

            _particles.parallelForEach([&](const int i) {
                VectorDr pos_i = _particles.positions[i];

                _particles.forEachNearby(pos_i, [&](const int j, const VectorDr & pos_j) {
                    VectorDr pos_ij = pos_i - pos_j;
                    VectorDr v_ij   = _velocities[i] - _velocities[j];

                    _deltaVelocity[i] += _particles.mass() * 2 * _viscosityCoeff * h * _c_s
                        / (_particles.densities[i] + _particles.densities[j]) * std::min(v_ij.dot(pos_ij), 0.)
                        / (pos_ij.dot(pos_ij) + 0.01 * h * h) * _particles.gradientKernel(pos_ij);
                });
            });

            _particles.parallelForEach([&](const int i) { _velocities[i] += _deltaVelocity[i] * dt; });
        }
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::applySurfaceTension(const real dt) {
        if (_surfaceTension) {
            _particles.parallelForEach([&](const int i) {
                VectorDr pos_i = _particles.positions[i];
                VectorDr delta = VectorDr::Zero();

                _particles.forEachNearby(pos_i, [&](const int j, const VectorDr & pos_j) {
                    VectorDr pos_ij = pos_i - pos_j;
                    delta += _particles.kernel(pos_ij) * pos_ij;
                });

                _velocities[i] -= _surfaceTension * delta * dt;
            });
        }
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::applyPressure(const real dt) {
        _particles.parallelForEach([&](const int i) {
            VectorDr delta = VectorDr::Zero();
            VectorDr pos_i = _particles.positions[i];
            _particles.forEachNearby(pos_i, [&](const int j, const VectorDr & pos_j) {
                delta -= _particles.mass()
                    * (_pressures[i] / _particles.densities[i] / _particles.densities[i]
                       + _pressures[j] / _particles.densities[j] / _particles.densities[j])
                    * _particles.gradientKernel(pos_j - pos_i);
            });

            _velocities[i] += delta * dt;
        });
    }

    template<int Dim> void WeakCompSphLiquid<Dim>::updateDensity(const real dt) {
        _particles.parallelForEach([&](const int i) {
            VectorDr pos_i = _particles.positions[i];
            double   delta = 0;
            _particles.forEachNearby(pos_i, [&](const int j, const VectorDr & pos_j) {
                delta += _velocities[j].dot(_particles.gradientKernel(pos_i - pos_j));
            });

            _particles.densities[i] += delta * _particles.mass() * dt;
        });
    }

    template class WeakCompSphLiquid<2>;
    template class WeakCompSphLiquid<3>;

} // namespace PhysX
