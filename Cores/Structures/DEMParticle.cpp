#include "DEMParticle.h"

#include "Structures/StaggeredGrid.h"

#include "Solvers/IterativeSolver.h"

#include "Utilities/Types.h"

#include <numbers>

#include <cmath>

namespace PhysX {
    template<int Dim>
    DEMParticle<Dim>::DEMParticle(
        const real radius, const real mass, const size_t cnt, const VectorDr & pos):
        SmoothedParticles<Dim>(radius, cnt, pos, mass, 2){
        Young = 1e9;
        Poisson = 0.3;
        contact_angle = 30. / 180. * std::numbers::pi;
        tan_fricangle = std::tan(0.5);
        K_norm = Young * radius;
        K_tang = K_norm * Poisson;
        volume_liquid_bridge = 4. / 3. * std::numbers::pi * std::pow(radius, 3.) * 0.01 * 0.01;
        d_rupture = (1.f + 0.5f * contact_angle) * (std::pow(volume_liquid_bridge, 1. / 3.) + 0.1 * std::pow(volume_liquid_bridge, 2. / 3.));
        c0 = 0.8;
        cmc = 1.;
        cmcp = 0.1;
        csat = 0.1;
        G = QuadraticBezierCoeff(c0, cmc, cmcp, csat);
        sr = cmcp;
        printf("drupture: %lf\n", d_rupture/radius);
    }

    template<int Dim>
    real DEMParticle<Dim>::_tan_fricangle(){return tan_fricangle;}

    template<int Dim>
    real DEMParticle<Dim>::_K_norm(){return K_norm;}

    template<int Dim>
    real DEMParticle<Dim>::_K_tang(){return K_tang;}

    template<int Dim>
    void DEMParticle<Dim>::setYoung(const real k){
        Young = k;
        K_norm = Young * _radius;
        K_tang = K_norm * Poisson;
    }

    template<int Dim>
    void DEMParticle<Dim>::setPoisson(const real k){
        Poisson = k;
        K_tang = K_norm * Poisson;
    }

    template<int Dim>
    void DEMParticle<Dim>::setfricangle(const real k){
        tan_fricangle = std::tan(k);
    }

    template<int Dim> 
    auto DEMParticle<Dim>::getForce(const int i, const int j) const -> VectorDr {
        VectorDr dij = positions[j] - positions[i];
        VectorDr vij = velocities[j] - velocities[i];
        return ComputeDemForces(dij, vij) + ComputeDemCapillaryForces(dij, vij);
    }

    template<int Dim>
    auto DEMParticle<Dim>::ComputeDemForces(const VectorDr & dij, const VectorDr & vij) const -> VectorDr {
        VectorDr f = VectorDr::Zero();
        real dist = dij.norm();
        real penetration_depth = 2 * _radius - dist;
        if (penetration_depth > 0.)
        {
            VectorDr n = VectorDr::Zero();
            if(dist <= 0.0001 * _radius)
                dist = 0.0001 * _radius;

            n = dij * (1 / dist);
            real dot_epslion = vij.dot(n);
            VectorDr vij_tangential = vij - dot_epslion * n;

            VectorDr normal_force = K_norm * penetration_depth * n;
            VectorDr shear_force = -K_tang * vij_tangential;

            real max_fs = normal_force.norm() * tan_fricangle;
            real shear_force_norm = shear_force.norm();

            if (shear_force_norm > max_fs){
                shear_force = shear_force * max_fs / shear_force_norm;
            }
            f = -normal_force - shear_force;
        }
        return f;
    }

    template<int Dim>
    auto DEMParticle<Dim>::ComputeDemCapillaryForces(const VectorDr & dij, const VectorDr & vij) const -> VectorDr {
        VectorDr f = VectorDr::Zero();
        real dist = dij.norm();
        real H = dist - 2 * _radius;
        if (H < d_rupture && H > 0.000001)
        {
            VectorDr n = VectorDr::Zero();
            n = dij * (1 / dist);
            real dot_epslion = vij.dot(n);
            VectorDr vij_normal = dot_epslion * n;
            VectorDr vij_tangential = vij - dot_epslion * n;

            // float coeff_c = csat + (1.f - sr) * (c0 - csat);
            real coeff_c = G.calculate(sr);

            // printf("cohesive=%.3f \n", coeff_c);

            real d = -H + std::sqrt(H * H + volume_liquid_bridge / (std::numbers::pi * _radius));
            real phi = std::sqrt(2. * H / _radius * (-1.f + sqrtf(1.f + volume_liquid_bridge / (std::numbers::pi * _radius * H * H))));
            real neck_curvature_pressure = -2. * std::numbers::pi * coeff_c * _radius * std::cos(contact_angle) / (1. + H / (2. * d));
            real surface_tension_force = -2. * std::numbers::pi * coeff_c * _radius * phi * std::sin(contact_angle);

            f = -n * (neck_curvature_pressure + surface_tension_force);
        }
        return f;
    }

    template<int Dim>
    auto DEMParticle<Dim>::getForceSum(const int i) const -> VectorDr {
        VectorDr f = VectorDr::Zero();
        forEachNearby(
            positions[i], [&](const int j, const VectorDr & nearbyPos) { 
                f += getForce(i, j); 
            });
        return f;
    }



    
    template class DEMParticle<2>;
    template class DEMParticle<3>;

} // namespace PhysX
