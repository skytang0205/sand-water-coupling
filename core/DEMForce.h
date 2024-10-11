#pragma once

#include "Particle.h"
#include "Common.h"

namespace Pivot {
    class DEMForce{
        public:
        DEMForce(double radius){
            _radius = radius;
            Young = 1e6;
            Poisson = 0.1;
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
            //G = QuadraticBezierCoeff(c0, cmc, cmcp, csat);
            surface_tensor_cof = 1.0;
            //printf("drupture: %lf\n", d_rupture/radius);
        }

        double _tan_fricangle(){return tan_fricangle;}

        double _K_norm(){return K_norm;}

        double _K_tang(){return K_tang;}

        double _Young() const{return Young;}

        double _Poission() const{return Poisson;}

        void setYoung(double k){
            Young = k;
            K_norm = Young * _radius;
            K_tang = K_norm * Poisson;
        }

        void setPoisson(double k){
            Poisson = k;
            K_tang = K_norm * Poisson;
        }

        void setfricangle(double k){
            tan_fricangle = std::tan(k);
        }


        Vector2d getForce(Particle &particlei, Particle &particlej){
            Vector2d dij = particlej.Position - particlei.Position;
            Vector2d vij = particlej.Velocity - particlei.Velocity;
            //double sr = (particlei.SaturateRate + particlej.SaturateRate) * .5;
            return ComputeDemForces(dij, vij);// + ComputeDemCapillaryForces(dij, vij, sr);
        }
        Vector2d ComputeDemForces(Vector2d & dij, Vector2d & vij){
            Vector2d f = Vector2d::Zero();
            double dist = dij.norm();
            double penetration_depth = 2 * _radius - dist;
            if (penetration_depth > 0.)
            {
                Vector2d n = Vector2d::Zero();
                if(dist <= 0.0001 * _radius)
                    return Vector2d::Zero();

                n = dij * (1 / dist);
                double dot_epslion = vij.dot(n);
                Vector2d vij_tangential = vij - dot_epslion * n;

                Vector2d normal_force = K_norm * penetration_depth * n;
                Vector2d shear_force = -K_tang * vij_tangential;

                double max_fs = normal_force.norm() * tan_fricangle;
                double shear_force_norm = shear_force.norm();

                if (shear_force_norm > max_fs){
                    shear_force = shear_force * max_fs / shear_force_norm;
                }
                f = -normal_force - shear_force;
            }
            return f;
        }

        // Vector2d ComputeDemCapillaryForces(Vector2d & dij, Vector2d & vij, double sr){
        //     Vector2d f = Vector2d::Zero();
        //     double dist = dij.norm();
        //     double H = dist - 2 * _radius;
        //     if (H < d_rupture && H > 0.000001)
        //     {
        //         Vector2d n = Vector2d::Zero();
        //         n = dij * (1 / dist);
        //         double dot_epslion = vij.dot(n);
        //         Vector2d vij_normal = dot_epslion * n;
        //         Vector2d vij_tangential = vij - dot_epslion * n;

        //         // float coeff_c = csat + (1.f - sr) * (c0 - csat);
        //         double coeff_c = G.calculate(sr) * surface_tensor_cof;

        //         // printf("cohesive=%.3f \n", coeff_c);

        //         double d = -H + std::sqrt(H * H + volume_liquid_bridge / (std::numbers::pi * _radius));
        //         double phi = std::sqrt(2. * H / _radius * (-1.f + sqrtf(1.f + volume_liquid_bridge / (std::numbers::pi * _radius * H * H))));
        //         double neck_curvature_pressure = -2. * std::numbers::pi * coeff_c * _radius * std::cos(contact_angle) / (1. + H / (2. * d));
        //         double surface_tension_force = -2. * std::numbers::pi * coeff_c * _radius * phi * std::sin(contact_angle);

        //         f = -n * (neck_curvature_pressure + surface_tension_force);
        //     }
        //     return f;
        // }

        Vector2d getForceSum(GridData<std::vector<Particle>> const &grid, Particle & particle){
            Vector2d f = Vector2d::Zero();
            Vector2i const lower = grid.GetGrid().Clamp(grid.GetGrid().CalcLower<1>(particle.Position));
            int range = int(_radius * grid.GetGrid().GetInvSpacing() + 2.) * 2;
			Vector2i const size = grid.GetGrid().GetSize();
			for(int i = std::max(0, lower.x()-range); i < std::min(size.x(), lower.x()+range); i++){
				for(int j = std::max(0, lower.y()-range); j < std::min(size.y(), lower.y()+range); j++){
                    for(auto nb : grid[Vector2i(i,j)]){
                        if((nb.Position-particle.Position).squaredNorm() < 4 * _radius * _radius)
                            f += getForce(particle, nb);
                    }
			    }
			}
            return f;
        }

        Vector2d getForceSum(std::vector<Particle> const &m_Particles, Particle & p0){
            Vector2d f = Vector2d::Zero();
			for (auto p1 : m_Particles) {
				if (((p0.Position - p1.Position).squaredNorm() < 4 * _radius * _radius)) {
					f = f + getForce(p0, p1);
				}
			}
            return f;
        }

        private:
        double _radius;
        double Young, Poisson, K_norm, K_tang, tan_fricangle;
        double contact_angle, volume_liquid_bridge, d_rupture;
        double c0, cmc, cmcp, csat, sr, surface_tensor_cof;
        //QuadraticBezierCoeff G;
    };
}