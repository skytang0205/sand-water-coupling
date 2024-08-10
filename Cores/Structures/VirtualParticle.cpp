#include "VirtualParticle.h"

#include "Structures/StaggeredGrid.h"

#include "Solvers/IterativeSolver.h"

#include <numbers>

#include <cmath>

namespace PhysX {
    template<int Dim>
    VirtualParticle<Dim>::VirtualParticle(
        const StaggeredGrid<Dim> & grid, const real radius, const size_t cnt, const VectorDr & pos):
        _grid(grid),
        _weight(_grid.nodeGrid(), 0), SmoothedParticles<Dim>(radius, cnt, pos), _alpha_0(-1e6) {}

    template<int Dim> void VirtualParticle<Dim>::computeDensities(const SmoothedParticles<Dim> & realParticles) {
        densities._data.resize(positions.size());
        parallelForEach([&](const int I) {
            const VectorDr p_I = positions[I];
            densities[I]       = 0;
            realParticles.forEachNearby(
                p_I, [&](const int j, const VectorDr & p_j) { densities[I] += kernel(p_I - p_j); });
            densities[I] *= realParticles.mass();
        });
    }

    template<int Dim> void VirtualParticle<Dim>::computeVolumes(const SmoothedParticles<Dim> & realParticles) {
        volumes._data.resize(positions.size());
        parallelForEach([&](const int I) {
            const VectorDr p_I = positions[I];
            volumes[I]         = 0;
            realParticles.forEachNearby(p_I, [&](const int j, const VectorDr & p_j) {
                volumes[I] += realParticles.volumes[j] * realParticles.volumes[j] * realParticles.kernel(p_I - p_j);
            });
        });
    }

    template<int Dim> void VirtualParticle<Dim>::computeInfo(const SmoothedParticles<Dim> & realParticles) {
        computeDensities(realParticles);
        computeVolumes(realParticles);
    }

    template<int Dim> void VirtualParticle<Dim>::generateParticles(const SmoothedParticles<Dim> & realParticles) {
        clear();

        _weight.parallelForEach([&](const VectorDi & I) {
            VectorDr p_I = _weight.position(I);
            _weight[I]   = 0;
            realParticles.forEachNearby(
                p_I, [&](int j, const VectorDr & p_j) { _weight[I] += realParticles.kernel(p_I - p_j); });
        });

        _weight.forEach([&](const VectorDi & I) {
            if (_weight[I] > 1e-6) add(_weight.position(I));
        });

        resetNearbySearcher();

        computeInfo(realParticles);
    }

    template<int Dim>
    void VirtualParticle<Dim>::setAlpha(const SmoothedParticles<Dim> & realParticles, const real target_rho) {
        generateParticles(realParticles);

        calculateCoef(target_rho);

        _alpha_0 = _matLaplacian.diagonal().maxCoeff();

        std::cout << _alpha_0 << std::endl;
    }

    template<int Dim>
    void VirtualParticle<Dim>::project(
        ParticlesBasedVectorField<Dim> &       velocity,
        const ParticlesBasedVectorField<Dim> & boundary_velocity,
        const SmoothedParticles<Dim> &         real_particle,
        const SmoothedParticles<Dim> &         boundary_particle,
        const real                             target_rho) {
        buildLinearSystem(velocity, boundary_velocity, real_particle, boundary_particle, target_rho);
        solveLinearSystem();
        applyPressureGradient(velocity, real_particle, target_rho);
    }

    template<int Dim>
    void VirtualParticle<Dim>::buildLinearSystem(
        const ParticlesBasedVectorField<Dim> & velocity,
        const ParticlesBasedVectorField<Dim> & boundary_velocity,
        const SmoothedParticles<Dim> &         real_particle,
        const SmoothedParticles<Dim> &         boundary_particle,
        const real                             target_rho) {
        _pressures.resize(this);
        _divergence.resize(this);
        _velocities.resize(this);

        calculateVelocity(velocity, boundary_velocity, real_particle, boundary_particle);
        calculateDiv(velocity, boundary_velocity, real_particle, boundary_particle, target_rho);

        calculateCoef(target_rho);
    }

    template<int Dim> void VirtualParticle<Dim>::solveLinearSystem() {
        IterativeSolver::solve(_matLaplacian, _pressures.asVectorXr(), _divergence.asVectorXr());
    }

    template<int Dim>
    void VirtualParticle<Dim>::applyPressureGradient(
        ParticlesBasedVectorField<Dim> & velocity,
        const SmoothedParticles<Dim> &   real_particle,
        const real                       target_rho) {
        real_particle.parallelForEach([&](const int i) {
            VectorDr p_i   = real_particle.positions[i];
            VectorDr force = VectorDr::Zero();
            forEachNearby(p_i, [&](const int J, const VectorDr & p_J) {
                VectorDr r_iJ = p_i - p_J;
                force += volumes[J] * _pressures[J] * gradientKernel(r_iJ);
            });

            _velocities[i] -= force / target_rho;
        });
    }

    template<int Dim>
    void VirtualParticle<Dim>::calculateVelocity(
        const ParticlesBasedVectorField<Dim> & velocity,
        const ParticlesBasedVectorField<Dim> & boundary_velocity,
        const SmoothedParticles<Dim> &         real_particle,
        const SmoothedParticles<Dim> &         boundary_particle) {
        parallelForEach([&](const int I) {
            VectorDr p_I   = positions[I];
            _velocities[I] = (velocity(p_I) + boundary_velocity(p_I))
                / (real_particle.getBiasWeight(p_I) + boundary_particle.getBiasWeight(p_I) + 1e-6);
        });
    }

    template<int Dim> void VirtualParticle<Dim>::calculateCoef(const real target_rho) {
        _matLaplacian.resize(positions.size(), positions.size());
        _coefficients.clear();
        forEach([&](const int I) {
            VectorDr p_I = positions[I];
            double   sum = 0;
            forEachNearby(p_I, [&](int J, const VectorDr p_J) {
                double r_ij     = (p_I - p_J).norm();
                double alpha_ij = 2. / target_rho * volumes[J] * firstDerivativeKernel(r_ij) / (r_ij + 1e-6);
                sum -= alpha_ij;
                _coefficients.push_back({ I, J, alpha_ij });
            });
            //_coefficients.push_back({ I, I, sum });
            _coefficients.push_back({ I, I, std::max(sum, _alpha_0) });
        });

        _matLaplacian.setFromTriplets(_coefficients.begin(), _coefficients.end());

        std::cout << _matLaplacian.coeff(0, 0) << std::endl;
        std::cout << _matLaplacian.coeff(0, 1) << std::endl;
    }

    template<int Dim>
    void VirtualParticle<Dim>::calculateDiv(
        const ParticlesBasedVectorField<Dim> & velocity,
        const ParticlesBasedVectorField<Dim> & boundary_velocity,
        const SmoothedParticles<Dim> &         real_particle,
        const SmoothedParticles<Dim> &         boundary_particle,
        const real                             target_rho) {
        parallelForEach([&](const int I) {
            VectorDr p_I = positions[I];
            _divergence[I] =
                (velocity.divergence(p_I) - _velocities[I].dot(real_particle.getBiasGradient(p_I))
                 + boundary_velocity.divergence(p_I) - _velocities[I].dot(boundary_particle.getBiasGradient(p_I))
                 + kappa * std::max(densities[I] - target_rho, 0.) / target_rho);

            //_particles.forEachNearby(p_I, [&](const int j, const VectorDr & p_j) {
            //    VectorDr r_Ij = p_I - p_j;
            //    _divergence[I] += _particles.volumes[j] * (_velocities[j] - _virtual_velocity[I]).dot(r_Ij)
            //        / r_Ij.norm() * _particles.firstDerivativeKernel(r_Ij.norm());
            //});

            //_boundary_particles.forEachNearby(p_I, [&](const int j, const VectorDr & p_j) {
            //    VectorDr r_Ij = p_I - p_j;
            //    _divergence[I] += _boundary_particles.volumes[j]
            //        * (_boundary_velocity[j] - _virtual_velocity[I]).dot(r_Ij) / r_Ij.norm()
            //        * _boundary_particles.firstDerivativeKernel(r_Ij.norm());
            //});

            //_divergence[I] -= _particles._mass / _targetDensity
            //    * std::max(_virtual_particles.densities[I] - _targetDensity, 0.) / _targetDensity;
            //_divergence[I] /= dt;
        });
    }

    template class VirtualParticle<2>;
    template class VirtualParticle<3>;

} // namespace PhysX
