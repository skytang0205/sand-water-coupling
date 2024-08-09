#include "SmthParticleHydrodLiquid.h"

#include "Utilities/Constants.h"

#include "Solvers/IterativeSolver.h"

namespace PhysX {

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::writeDescription(YAML::Node & root) const {
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
            _particles.forEach([&](const int i) {
                IO::writeValue(fout, _particles.positions[i].template cast<float>().eval());
            });
            if constexpr (Dim == 3) {
                _particles.forEach([&](const int i) {
                    IO::writeValue(fout, VectorDf::Unit(2).eval());
                });
            }
            _particles.forEach([&](const int i) {
                IO::writeValue(fout, float(_velocities[i].norm()));
            });
        }
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::saveFrame(const std::string & frameDir) const {
        { // Save particles.
            std::ofstream fout(frameDir + "/particles.sav", std::ios::binary);
            _particles.positions.save(fout);
        }
        { // save velocities.
            std::ofstream fout(frameDir + "/velocities.sav", std::ios::binary);
            _velocities.save(fout);
        }
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::loadFrame(const std::string & frameDir) {
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

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::initialize() {
        _particles.resetNearbySearcher();

        _particles.computeDensities();
        _particles.computeVolumes();

        generateVirtualParticle();

        calculateCoef();

        _alpha_0 = _matLaplacian.diagonal().maxCoeff();

        double omega = 2.;

        if constexpr (Dim == 2)
            _particles.parallelForEach([&](const int i) {
                _velocities[i] = omega * _particles.positions[i].y() * VectorDr::Unit(0) - omega * _particles.positions[i].x() * VectorDr::Unit(1);
            });
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::advance(const real dt) {
        moveParticles(dt);

        applyExternalForces(dt);
        applyViscosityForce(dt);

        generateVirtualParticle();

        calculateVirtualVelocity();

        applyPressureForce(dt);
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::reinitializeParticlesBasedData() {
        _pressures.resize(&_virtual_particles);
        _divergence.resize(&_virtual_particles);
        _virtual_velocity.resize(&_virtual_particles);
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::moveParticles(const real dt) {
        _particles.parallelForEach([&](const int i) {
            _particles.positions[i] += _velocities[i] * dt;
        });

        // Resolve collisions.
        for (const auto & collider : _colliders) {
            collider->collide(_particles.positions, _velocities, _particles.radius());
        }

        _particles.resetNearbySearcher();
        _particles.computeDensities();
        _particles.computeVolumes();
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::applyExternalForces(const real dt) {
        if (_enableGravity) {
            _particles.parallelForEach([&](const int i) {
                _velocities[i][1] -= kGravity * dt;
            });
        }
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::applyViscosityForce(const real dt) {
        if (_viscosityCoeff) {
            auto newVelocities = _velocities;
            _particles.parallelForEach([&](const int i) {
                newVelocities[i] += _viscosityCoeff * _velocities.divFreeLaplacianAtDataPoint(i) * dt;
            });
            _velocities = newVelocities;
        }
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::applyPressureForce(const real dt) {
        //// Compute pressures by the equation of state.
        //_particles.parallelForEach([&](const int i) {
        //    _pressures[i] = _eosMultiplier * (_particles.densities[i] - _targetDensity);
        //    if (_pressures[i] < 0) _pressures[i] = 0;
        //});

        //// Apply pressure gradient.
        //_particles.parallelForEach([&](const int i) {
        //    _velocities[i] -= _pressures.symmetricGradientAtDataPoint(i) / _particles.densities[i] * dt;
        //});
        calculateCoef();
        calculateDiv(dt);
        IterativeSolver::solve(
            _matLaplacian,
            _pressures.asVectorXr(),
            _divergence.asVectorXr());

        _particles.parallelForEach([&](const int i) {
            VectorDr p_i   = _particles.positions[i];
            VectorDr force = VectorDr::Zero();
            _virtual_particles.forEachNearby(p_i, [&](const int J, const VectorDr & p_J) {
                VectorDr r_iJ = p_i - p_J;
                force += _virtual_particles.volumes[J] * _pressures[J] * r_iJ / r_iJ.norm() * _virtual_particles.firstDerivativeKernel(r_iJ.norm());
            });
            _velocities[i] -= force * dt / _targetDensity;
        });
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::generateVirtualParticle() {
        _virtual_weight.parallelForEach([&](const VectorDi I) {
            VectorDr p_I       = _virtual_weight.position(I);
            _virtual_weight[I] = 0;
            _particles.forEachNearby(p_I, [&](int j, const VectorDr p_j) {
                _virtual_weight[I] += _particles.kernel(p_I - p_j);
            });
        });

        _virtual_particles.positions._data.clear();
        _virtual_weight.forEach([&](const VectorDi I) {
            if (_virtual_weight[I] > 1e-6)
                _virtual_particles.positions._data.push_back(_virtual_weight.position(I));
        });

        _virtual_particles.resetNearbySearcher();
        reinitializeParticlesBasedData();
        _virtual_particles.computeVirtualVolumes(_particles);
        _virtual_particles.computeVirtualDensities(_particles);
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::calculateCoef() {
        _matLaplacian.resize(_virtual_particles.size(), _virtual_particles.size());
        _coefficients.clear();
        _virtual_particles.forEach([&](const int I) {
            VectorDr p_I = _virtual_particles.positions[I];
            double   sum = 0;
            _virtual_particles.forEachNearby(p_I, [&](int J, const VectorDr p_J) {
                double r_ij     = (p_I - p_J).norm();
                double alpha_ij = 2. / _targetDensity * _virtual_particles.volumes[J] * _virtual_particles.firstDerivativeKernel(r_ij) / (r_ij + 1e-6);
                sum += alpha_ij;
                _coefficients.push_back({ I, J, alpha_ij });
            });
            _coefficients.push_back({ I, I, std::max(sum, -_alpha_0) });
        });

        _matLaplacian.setFromTriplets(_coefficients.begin(), _coefficients.end());
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::calculateVirtualVelocity() {
        _virtual_particles.parallelForEach([&](const int I) {
            VectorDr p_I = _virtual_particles.positions[I];
            double   sum = 0;
            _virtual_velocity[I].setZero();
            _particles.forEachNearby(p_I, [&](const int j, const VectorDr & p_j) {
                _virtual_velocity[I] += _particles.volumes[j] * _velocities[j] * _particles.kernel(p_I - p_j);
                sum += _particles.volumes[j] * _particles.kernel(p_I - p_j);
            });
            _virtual_velocity[I] /= (sum + 1e-6);
        });
    }

    template<int Dim>
    void SmthParticleHydrodLiquid<Dim>::calculateDiv(real dt) {
        _virtual_particles.parallelForEach([&](const int I) {
            VectorDr p_I   = _virtual_particles.positions[I];
            _divergence[I] = 0;
            _particles.forEachNearby(p_I, [&](const int j, const VectorDr & p_j) {
                VectorDr r_Ij = _virtual_particles.positions[I] - _particles.positions[j];
                _divergence[I] += _particles.volumes[j] * (_velocities[j] - _virtual_velocity[I]).dot(r_Ij) / r_Ij.norm() * _particles.firstDerivativeKernel(r_Ij.norm());
            });
            _divergence[I] -= _particles._mass / _targetDensity * std::max(_virtual_particles.densities[I] - _targetDensity, 0.) / _targetDensity;
            _divergence[I] /= dt;
        });
    }

    template class SmthParticleHydrodLiquid<2>;
    template class SmthParticleHydrodLiquid<3>;

} // namespace PhysX
