#pragma once

#include "Geometries/ImplicitSurface.h"
#include "Physics/SpringMassSystem.h"

#include <fmt/core.h>

#include <numbers>

namespace PhysX {

    class SpringMassSystemBuilder final {
    public:
        template<int Dim>
        static std::unique_ptr<SpringMassSystem<Dim>> build(const int option) {
            switch (option) {
            case 0:
                return buildCase0<Dim>();
            case 1:
                return buildCase1<Dim>();
            case 2:
                return buildCase2<Dim>();
            default:
                reportError("invalid option");
                return nullptr;
            }
        }

    protected:
        template<int Dim>
        static std::unique_ptr<SpringMassSystem<Dim>> buildCase0() {
            DECLARE_DIM_TYPES(Dim)
            auto smSystem            = std::make_unique<SpringMassSystem<Dim>>();
            smSystem->_enableGravity = false;

            smSystem->_particles.add(VectorDr::Zero());
            smSystem->_particles.add(VectorDr::Unit(0));

            smSystem->_velocities.resize(&smSystem->_particles);

            smSystem->_springs.push_back({ 0, 1, real(.5), 100, 1 });

            for (int axis = 0; axis < Dim; axis++)
                smSystem->_constrainedDofs.insert(axis);

            return smSystem;
        }

        template<int Dim>
        static std::unique_ptr<SpringMassSystem<Dim>> buildCase1() {
            DECLARE_DIM_TYPES(Dim)
            auto smSystem = std::make_unique<SpringMassSystem<Dim>>();

            for (int i = 0; i < 10; i++)
                smSystem->_particles.add(VectorDr::Unit(0) * i / 4);

            smSystem->_velocities.resize(&smSystem->_particles);

            for (int i = 0; i < 9; i++)
                smSystem->_springs.push_back({ i, i + 1, real(.25), 1000, 0 });

            smSystem->_colliders.push_back(
                std::make_unique<StaticCollider<Dim>>(
                    std::make_unique<ImplicitPlane<Dim>>(VectorDr::Zero() - VectorDr::Unit(1) * 2, VectorDr::Unit(1)),
                    real(.5)));

            for (int axis = 0; axis < Dim; axis++)
                smSystem->_constrainedDofs.insert(axis);

            return smSystem;
        }

        template<int Dim>
        static std::unique_ptr<SpringMassSystem<Dim>> buildCase2() {
            DECLARE_DIM_TYPES(Dim)
            auto smSystem = std::make_unique<SpringMassSystem<Dim>>();

            const int  xScale         = 10;
            const int  yScale         = 2;
            const real length         = real(2.5);
            const real spacing        = length / xScale;
            const real stiffnessCoeff = 10000;
            const real dampingCoeff   = 0.1;

            const auto id = [=](const int i, const int j) { return j * (xScale + 1) + i; };

            for (int j = 0; j <= yScale; j++)
                for (int i = 0; i <= xScale; i++)
                    smSystem->_particles.add((VectorDr::Unit(0) * i + VectorDr::Unit(1) * j) * spacing);

            smSystem->_velocities.resize(&smSystem->_particles);

            for (int j = 0; j <= yScale; j++)
                for (int i = 0; i < xScale; i++)
                    smSystem->_springs.push_back({ id(i, j), id(i + 1, j), spacing, stiffnessCoeff, dampingCoeff });

            for (int j = 0; j < yScale; j++)
                for (int i = 0; i <= xScale; i++)
                    smSystem->_springs.push_back({ id(i, j), id(i, j + 1), spacing, stiffnessCoeff, dampingCoeff });

            for (int j = 0; j < yScale; j++)
                for (int i = 0; i < xScale; i++)
                    smSystem->_springs.push_back({ id(i, j), id(i + 1, j + 1), spacing * real(std::numbers::sqrt2), stiffnessCoeff, dampingCoeff });

            for (int j = 1; j <= yScale; j++)
                for (int i = 0; i < xScale; i++)
                    smSystem->_springs.push_back({ id(i, j), id(i + 1, j - 1), spacing * real(std::numbers::sqrt2), stiffnessCoeff, dampingCoeff });

            for (int j = 0; j <= yScale; j++)
                for (int axis = 0; axis < Dim; axis++)
                    smSystem->_constrainedDofs.insert(id(0, j) * Dim + axis);

            return smSystem;
        }

        static void reportError(const std::string & msg) {
            std::cerr << "Error: [SpringMassSystemBuilder] encountered " << msg << ".\n"
                      << msg << std::endl;

            std::exit(-1);
        }
    };

} // namespace PhysX
