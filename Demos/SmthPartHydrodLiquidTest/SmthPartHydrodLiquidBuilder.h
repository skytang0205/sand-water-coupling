#pragma once

#include "Geometries/ImplicitSurface.h"
#include "Physics/PredCorrIncomprSphLiquid.h"
#include "Physics/SmthParticleHydrodLiquid.h"
#include "Structures/StaggeredGrid.h"

#include <fmt/core.h>

#include <memory>

namespace PhysX {

    class SmthPartHydrodLiquidBuilder final {
    public:
        template<int Dim>
        static std::unique_ptr<SmthParticleHydrodLiquid<Dim>> build(const int scale, const int option, const bool pci) {
            switch (option) {
            case 0: return buildCase0<Dim>(scale, pci);
            default: reportError("invalid option"); return nullptr;
            }
        }

    protected:
        template<int Dim> static std::unique_ptr<SmthParticleHydrodLiquid<Dim>> buildCase0(int scale, const bool pci) {
            DECLARE_DIM_TYPES(Dim)
            if (scale < 0) scale = 30;
            const real length = real(2);

            const VectorDi     resolution = 3 * scale * VectorDi::Ones();
            StaggeredGrid<Dim> grid(2, length / scale / 2, resolution);
            const real         density = 1000;
            const real         radius  = length / 2 / scale / 2;
            auto               liquid  = makeLiquid<Dim>(grid, radius, pci);
            liquid->_particles.generateBoxPacked(VectorDr::Zero(), VectorDr::Ones() * length / 4);
            liquid->_particles.setMass(density / liquid->_particles.getPackedKernelSum());
            liquid->_targetDensity = density;

            liquid->_velocities.resize(&liquid->_particles);

            liquid->_virtual_particles.setKappa(1. / liquid->_particles.getPackedKernelSum());

            // liquid->_colliders.push_back(
            //     std::make_unique<StaticCollider<Dim>>(
            //         std::make_unique<ComplementarySurface<Dim>>(
            //             std::make_unique<ImplicitBox<Dim>>(-length / 2 * VectorDr::Ones(), length *
            //             VectorDr::Ones()))));
            liquid->_boundary_particles.addSurface(
                ComplementarySurface<Dim>(
                    std::make_unique<ImplicitBox<Dim>>(-length / 2 * VectorDr::Ones(), length * VectorDr::Ones())),
                VectorDr::Zero(),
                VectorDr::Ones() * length * 1.2);

            liquid->_boundary_velocity.resize(&liquid->_boundary_particles);

            return liquid;
        }

        template<int Dim>
        static std::unique_ptr<SmthParticleHydrodLiquid<Dim>>
            makeLiquid(const StaggeredGrid<Dim> & grid, const real radius, const bool pci) {
            if (pci) return std::make_unique<PredCorrIncomprSphLiquid<Dim>>(grid, radius);
            else return std::make_unique<SmthParticleHydrodLiquid<Dim>>(grid, radius);
        }

        static void reportError(const std::string & msg) {
            std::cerr << "Error: [SmthPartHydrodLiquidBuilder] encountered " << msg << ".\n" << msg << std::endl;

            std::exit(-1);
        }
    };

} // namespace PhysX
