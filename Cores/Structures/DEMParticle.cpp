#include "DEMParticle.h"

#include "Structures/StaggeredGrid.h"

#include "Solvers/IterativeSolver.h"

#include <numbers>

#include <cmath>

namespace PhysX {
    template<int Dim>
    DEMParticle<Dim>::DEMParticle(
        const real radius, const real mass, const size_t cnt, const VectorDr & pos):
        SmoothedParticles<Dim>(radius, cnt, pos, mass, 4){
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

    
    template class DEMParticle<2>;
    template class DEMParticle<3>;

} // namespace PhysX
