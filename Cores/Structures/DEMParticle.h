#pragma once

#include "Structures/GridBasedScalarField.h"
#include "Structures/Particles.h"
#include "Structures/ParticlesNearbySearcher.h"
#include "Structures/SmoothedParticles.h"
#include "Structures/StaggeredGrid.h"

#include <iostream>
#include <cmath>

namespace PhysX {

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
        real K_norm    = 10000000;
        real K_tang    = 1000;
        real fricangle = 0.1;

    public:
        DEMParticle(
            const real radius, const size_t cnt = 0, const VectorDr & pos = VectorDr::Zero(), const real mass = 1):
            SmoothedParticles<Dim>(radius, cnt, pos, mass, 2){}


        DEMParticle & operator=(const DEMParticle & rhs) = delete;
        virtual ~DEMParticle()                            = default;

        void setK_norm(const real k) { K_norm = k; }
        void setK_tang(const real k) { K_tang = k; }
        void setfricangleg(const real k) { fricangle = k; }

        VectorDr getForce(const int i, const int j) const {
            VectorDr dis_pos = positions[i] - positions[j];
            real dis_pos_len = dis_pos.norm();
            VectorDr dis_pos_norm = VectorDr::Zero();
            if(dis_pos_len >=0.000001)
                dis_pos_norm = dis_pos / dis_pos_len;
            VectorDr dis_vel = velocities[i] - velocities[j];
            VectorDr F_norm = K_norm * (2 * _radius - dis_pos_len) * dis_pos_norm;
            VectorDr F_tang = - K_tang * (dis_vel - dis_vel.dot(dis_pos_norm) * dis_pos_norm);
            if(F_tang.norm() <= 0.000001)
                return F_norm + F_tang;
            return F_norm + (F_norm.norm()*fricangle > F_tang.norm() ?  F_tang : (F_norm.norm() * fricangle / F_tang.norm()) * F_tang);
        }

        VectorDr getForceSum(const int i) const {
            VectorDr f = VectorDr::Zero();
            forEachNearby(
                positions[i], [&](const int j, const VectorDr & nearbyPos) { 
                    f += getForce(i, j); 
                });
            return f;
        }

        void applyForce(
            const ParticlesBasedVectorField<Dim> & boundary_velocity,
            const BoundaryParticles<Dim> &         boundary_particle,
            const real                             dt) {
            parallelForEach([&](const int i) {
                VectorDr p_i   = positions[i];
                VectorDr force = VectorDr::Zero();
                force += getForceSum(i) * dt;
                //std::cout << getForceSum(i) << std::endl;

                //force = VectorDr::Zero(); 

                boundary_particle.forEachNearby(p_i, [&](const int js, const VectorDr & p_js) {
                    VectorDr dis_vel = velocities[i] - boundary_velocity[js];
                    real dis_vel_norm = boundary_particle.norms[js].dot(dis_vel);
                    VectorDr F_norm = - boundary_particle.volumes[js]
                        * std::min(dis_vel_norm, 0.)
                        * boundary_particle.norms[js] * boundary_particle.kernel(p_i - p_js);
                    VectorDr F_tang = - K_tang * (dis_vel - dis_vel_norm * boundary_particle.norms[js]);
                    if(F_tang.norm() <= 0.000001)
                        force += F_norm + F_tang;
                    else
                        force += F_norm + (F_norm.norm()*fricangle > F_tang.norm() ?  F_tang : (F_norm.norm() * fricangle / F_tang.norm()) * F_tang);
                  });

                velocities[i] += force;
            });
        }

    };

} // namespace PhysX
