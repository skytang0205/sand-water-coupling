#pragma once

#include "Geometries/Surface.h"
#include "Structures/Particles.h"
#include "Structures/ParticlesNearbySearcher.h"
#include "Structures/GridBasedScalarField.h"
#include "Structures/SmoothedParticles.h"
#include "Structures/StaggeredGrid.h"
#include "Structures/ParticlesBasedData.h"
#include "Utilities/Types.h"

namespace PhysX {

    class QuadraticBezierCoeff{   
    private:
        real a, b, c, d, e;
        real py0, py1, px1, py2;
    public:
       
        QuadraticBezierCoeff(const real py0 = 0., const real py1 = 0., const real px1 = 0., const real py2 = 0.): 
            py0(py0), py1(py1), px1(px1), py2(py2){
            a = (px1 + 1.f) / 4.f;
            b = -2.f * a;
            c = b * b;
            d = -4.f * (1.f + b - px1);
            e = 2.f * (1.f + b - px1);
        }

        real rx2t(const real sr) const
        {
            return (b + std::sqrt(c + d * (px1 - sr))) / e;
        }

        real calculate(const real sr) const
        {
            if (sr < 0.f)
                return py0;

            if (sr >= 1.f)
                return py2;

            if (sr <= px1)
            {
                const real t = sr / px1;
                const real omt = 1. - t;
                return omt * omt * py0 + 2 * t * omt * py1 + t * t * py1;
            }
            else
            {
                const real t = rx2t(sr);
                const real omt = 1. - t;
                return omt * omt * py1 + 2 * t * omt * py1 + t * t * py2;
            }
        }
    };

    template<int Dim> class DEMParticle : public SmoothedParticles<Dim> {
        DECLARE_DIM_TYPES(Dim)

    public:
        using SmoothedParticles<Dim>::positions;
        ParticlesBasedVectorData<Dim> velocities;
    protected:
        using SmoothedParticles<Dim>::_radius;
        using SmoothedParticles<Dim>::_kernelRadius;
        using SmoothedParticles<Dim>::_nearbySearcher;
        

    private:
        using Particles<Dim>::_mass;
        using Particles<Dim>::_invMass;

    public:
        using SmoothedParticles<Dim>::resetNearbySearcher;
        using SmoothedParticles<Dim>::setMass;
        using SmoothedParticles<Dim>::mass;
        using SmoothedParticles<Dim>::parallelForEach;
        using SmoothedParticles<Dim>::forEachNearby;
        using SmoothedParticles<Dim>::forEach;
        using SmoothedParticles<Dim>::radius;
        using SmoothedParticles<Dim>::kernelRadius;
        //using SmoothedParticles<Dim>::firstDerivativeKernel;

    private:
        real Young, Poisson, K_norm, K_tang, tan_fricangle;
        real contact_angle, volume_liquid_bridge, d_rupture;
        real c0, cmc, cmcp, csat, sr;
        QuadraticBezierCoeff G;

    public:
        DEMParticle(
            const real radius, const real mass = 1, const size_t cnt = 0, const VectorDr & pos = VectorDr::Zero());
        DEMParticle & operator=(const DEMParticle & rhs) = delete;
        virtual ~DEMParticle()                           = default;

        real _tan_fricangle();
        real _K_norm();
        real _K_tang();

        void setYoung(const real k);
        void setPoisson(const real k);
        void setfricangle(const real k);


        VectorDr getForce(const int i, const int j) const {
            VectorDr dij = positions[j] - positions[i];
            VectorDr vij = velocities[j] - velocities[i];
            return ComputeDemForces(dij, vij) ;//+ ComputeDemCapillaryForces(dij, vij);
        }

        VectorDr ComputeDemForces(const VectorDr & dij, const VectorDr & vij) const {
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
        
        VectorDr ComputeDemCapillaryForces(const VectorDr & dij, const VectorDr & vij) const {
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


        VectorDr getForceSum(const int i) const{
            VectorDr f = VectorDr::Zero();
            forEachNearby(
                positions[i], [&](const int j, const VectorDr & nearbyPos) { 
                    f += getForce(i, j); 
                });
            return f;
        }


    };

} // namespace PhysX
