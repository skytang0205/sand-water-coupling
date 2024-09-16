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
        real c0, cmc, cmcp, csat, sr, surface_tensor_cof;
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


        VectorDr getForce(const int i, const int j) const;

        VectorDr ComputeDemForces(const VectorDr & dij, const VectorDr & vij) const;
        
        VectorDr ComputeDemCapillaryForces(const VectorDr & dij, const VectorDr & vij) const;

        VectorDr getForceSum(const int i) const;
    };

} // namespace PhysX
