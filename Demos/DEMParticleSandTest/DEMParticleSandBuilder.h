#pragma once

#include "Geometries/ImplicitSurface.h"
#include "Physics/DEMParticleSand.h"
#include "Structures/StaggeredGrid.h"
#include "Utilities/Shapes.h"

#include <fmt/core.h>

#include <memory>

namespace PhysX {

    class DEMParticleSandBuilder final {
    public:
        template<int Dim>
        static std::unique_ptr<DEMParticleSand<Dim>> build(const int scale, const int option) {
            switch (option) {
            case 0: return buildCase0<Dim>(scale);
            default: reportError("invalid option"); return nullptr;
            }
        }

    protected:
        template<int Dim> static std::unique_ptr<DEMParticleSand<Dim>> buildCase0(int scale) {
            DECLARE_DIM_TYPES(Dim)
            if (scale < 0) scale = 30;
            const real length = real(1);

            const real         density = 2700;
            const real         radius  = length / 2 / scale / 2;
            auto               sand  = std::make_unique<DEMParticleSand<Dim>>(radius);
            auto               shape   = Shapes<Dim>(radius);           
            const real         omega   = 2.;
            const bool      if_Possion = false;
            VectorDr           box     = VectorDr::Ones() * length / 4;
            box(0) /= 2;
            shape.generateBox(VectorDr::Zero(), box / 2 , if_Possion);
            //shape.generateRotate(omega);
            sand->addShape(shape);
            //liquid->_particles.generateBoxPacked(VectorDr::Zero(), VectorDr::Ones() * length / 4);
            sand->_particles.setMass(density / sand->_particles.positions.size());

            // liquid->_colliders.push_back(
            //     std::make_unique<StaticCollider<Dim>>(
            //         std::make_unique<ComplementarySurface<Dim>>(
            //             std::make_unique<ImplicitBox<Dim>>(-length / 2 * VectorDr::Ones(), length *
            //             VectorDr::Ones()))));
            sand->_boundary_particles.addSurface(
                ComplementarySurface<Dim>(
                    std::make_unique<ImplicitBox<Dim>>(-length / 2 * VectorDr::Ones(), length * VectorDr::Ones())),
                VectorDr::Zero(),
                VectorDr::Ones() * length * 1.2);

            sand->_boundary_particles.setMass(density / sand->_particles.positions.size());

            sand->_boundary_velocity.resize(&sand->_boundary_particles);

            return sand;
        }


        static void reportError(const std::string & msg) {
            std::cerr << "Error: [DEMParticleSandBuilder] encountered " << msg << ".\n" << msg << std::endl;

            std::exit(-1);
        }
    };

} // namespace PhysX
