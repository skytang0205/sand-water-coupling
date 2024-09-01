#pragma once

#include "Structures/GridBasedScalarField.h"
#include "Structures/Particles.h"
#include "Structures/ParticlesNearbySearcher.h"
#include "Structures/SmoothedParticles.h"
#include "Structures/StaggeredGrid.h"

#include <cmath>

namespace PhysX {

    template<int Dim> class DEMParticle : public SmoothedParticles<Dim> {
        DECLARE_DIM_TYPES(Dim)

    public:
        using SmoothedParticles<Dim>::positions;
        ParticlesBasedVectorData<Dim>::velosities;
    protected:
        using SmoothedParticles<Dim>::_radius;
        using SmoothedParticles<Dim>::_kernelRadius;
        using SmoothedParticles<Dim>::_nearbySearcher;
        

    private:
        using SmoothedParticles<Dim>::_mass;
        using SmoothedParticles<Dim>::_invMass;

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
        const StaggeredGrid<Dim> _grid;
        GridBasedScalarData<Dim> _weight;

        real K_norm;
        real K_tang;
        real fricangle;

    public:
        DEMParticles(
            const real radius, const size_t cnt = 0, const VectorDr & pos = VectorDr::Zero(), const real mass = 1):
            Particles<Dim>(cnt, pos, mass),
            _radius(radius),_kernelRadius(radius * 2),
            _nearbySearcher(std::make_unique<HashGridSearcher<Dim>>(_kernelRadius)) {}

        DEMParticles & operator=(const DEMParticles & rhs) = delete;
        virtual ~DEMParticles()                            = default;

        void setK_norm(const real k) { K_norm = k; }
        void setK_tang(const real k) { K_tang = k; }
        void setfricangleg(const real k) { fricangle = k; }

        VectorDr getForce(const int i, const int j) const {
            VectorDr dis_pos = positions[i] - positions[j];
            real dis_pos_len = dis_pos.norm();
            VectorDr dis_pos_norm = dis_pos / dis_pos_len;
            VectorDr dis_vel = velosities[i] - velosities[j];
            VectorDr F_norm = K_norm * (2 * _radius - dis_pos_len) * dis_pos_norm;
            VectorDr F_tang = - K_tang * (dis_vel - dis_vec.dot(dis_pos_norm) * dis_pos_norm);
            return F_norm.norm()*fricangle > F_tang.norm() ? F_norm + F_tang : F_norm + F_norm.norm() * fricangle / F_tang.norm() * F_tang;
        }

        VectorDr getForceSum(const int i) const {
            VectorDr f = VectorDr::Zero();
            forEachNearby(
                positions[i], [&](const int j, const VectorDr & nearbyPos) { f += getForce(i, j); });
            return f;
        }

        void applyForce(
            const ParticlesBasedVectorField<Dim> & boundary_velocity,
            const BoundaryParticles<Dim> &         boundary_particle,
            const real                             dt) {
            parallelForEach([&](const int i) {
                VectorDr p_i   = positions[i];
                VectorDr force = getForceSum(i) * _invMass * dt;

                boundary_particle.forEachNearby(p_i, [&](const int js, const VectorDr & p_js) {
                    force -= boundary_particle.volumes[js]
                        * std::min(boundary_particle.norms[js].dot(velocity[i] - boundary_velocity[js]), 0.)
                        * boundary_particle.norms[js] * boundary_particle.kernel(p_i - p_js);
                });

                velocities[i] += force;
            });
        }

    };

} // namespace PhysX
