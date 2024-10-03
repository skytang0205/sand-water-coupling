#pragma once

#include "Structures/ParticlesAttribute.h"
#include "Structures/ParticlesBasedData.h"
#include "Structures/StaggeredGrid.h"
#include "Utilities/Types.h"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <vector>

namespace PhysX {
    template<int Dim> class PoissonDiskSampler {
        DECLARE_DIM_TYPES(Dim)
    public:
        PoissonDiskSampler(VectorDr boundary, VectorDr center, real radius):
            boundary(boundary), center(center), up_bound(center + boundary), low_bound(center - boundary),
            radius(radius) {
            gridCellSize = radius / std::sqrt(Dim);
            grid_size    = ((up_bound - low_bound) / gridCellSize).cast<int>();
            int size     = 1;
            for (int i = 0; i < Dim; i++) size *= grid_size(i);
            grid.resize(size, -1);
        }

        VectorDr generateNewPoint(VectorDr point) {
            static std::random_device               rd;
            static std::mt19937                     gen(rd());
            static std::uniform_real_distribution<> disX(low_bound(0), up_bound(0));
            static std::uniform_real_distribution<> disY(low_bound(1), up_bound(1));
            static std::uniform_real_distribution<> dis(0.0, 1.0);

            if constexpr (Dim == 2) {
                if (point == VectorDr::Zero()) {
                    VectorDr initialPoint = VectorDr::Ones();
                    initialPoint(0)       = disX(gen);
                    initialPoint(1)       = disY(gen);
                    return initialPoint;
                }
                real     angle    = dis(gen) * 2 * std::numbers::pi;
                real     dist     = radius + dis(gen) * radius;
                VectorDr newPoint = { std::cos(angle), std::sin(angle) };
                newPoint          = newPoint * dist + point;
                return newPoint;
            } else if constexpr (Dim == 3) {
                static std::uniform_real_distribution<> disZ(low_bound(2), up_bound(2));

                if (point == VectorDr::Zero()) {
                    VectorDr initialPoint = VectorDr::Ones();
                    initialPoint(0)       = disX(gen);
                    initialPoint(1)       = disY(gen);
                    initialPoint(2)       = disZ(gen);
                    return initialPoint;
                }
                real            angle1   = disX(gen) * 2 * std::numbers::pi;
                real            angle2   = std::acos(1 - 2 * disX(gen));
                real            dist     = radius + disX(gen) * radius * 0.5;
                Vector<3, real> newPoint = { std::cos(angle1) * std::sin(angle2),
                                             std::sin(angle1) * std::sin(angle2),
                                             std::cos(angle2) };

                newPoint = newPoint * dist + point;
                return newPoint;
            }
        }

        std::vector<VectorDr> generatePoints() {
            std::vector<VectorDr> samples;
            std::vector<VectorDr> processList;
            VectorDr              initialPoint = generateNewPoint(VectorDr::Zero());
            processList.push_back(initialPoint);
            addToGrid(initialPoint);

            while (! processList.empty()) {
                int      idx   = rand() % processList.size();
                VectorDr point = processList[idx];
                bool     found = false;

                for (int i = 0; i < 300; ++i) {
                    VectorDr newPoint = generateNewPoint(point);
                    if (isInBounds(newPoint) && isFarEnough(newPoint)) {
                        // samples.push_back(newPoint);
                        processList.push_back(newPoint);
                        addToGrid(newPoint);
                        found = true;
                        break;
                    }
                    // if(!isInBounds(newPoint)) printf("not in bounds\n");
                    // if (!isFarEnough(newPoint)) printf("not far\n");
                }
                if (! found) {
                    processList.erase(processList.begin() + idx);
                    // printf("yes\n");
                }
                printf("%zd %zd %zd\n", grid.size(), samples.size(), processList.size());
                if (gridPoints.size() > grid.size() * 1.1) break;
            }

            for (const auto & idx : grid) {
                if (idx != -1) { samples.push_back(gridPoints[idx]); }
            }
            return samples;
        }

        bool isInBounds(const VectorDr & point) const {
            VectorDi grid_idx = ((point - low_bound) / gridCellSize).cast<int>();
            if constexpr (Dim == 2)
                return grid_idx.x() >= 0 && grid_idx.x() < grid_size.x() && grid_idx.y() >= 0
                    && grid_idx.y() < grid_size.y();
            else if constexpr (Dim == 3)
                return grid_idx.x() >= 0 && grid_idx.x() < grid_size.x() && grid_idx.y() >= 0
                    && grid_idx.y() < grid_size.y() && grid_idx.z() >= 0 && grid_idx.z() < grid_size.z();
        }

    private:
        real                  radius, gridCellSize;
        VectorDr              boundary, center, low_bound, up_bound;
        VectorDi              grid_size;
        std::vector<int>      grid;
        std::vector<VectorDr> gridPoints;

        bool isFarEnough(const VectorDr & point) const {
            VectorDi grid_idx = ((point - low_bound) / gridCellSize).cast<int>();

            // for (const auto & p : gridPoints)
            //     if ((p - point).norm() < radius)
            //         return false;
            // return true;

            if constexpr (Dim == 2) {
                for (int y = std::max(0, grid_idx.y() - 3); y <= std::min(grid_size.y() - 1, grid_idx.y() + 3); ++y) {
                    for (int x = std::max(0, grid_idx.x() - 3); x <= std::min(grid_size.x() - 1, grid_idx.x() + 3);
                         ++x) {
                        int i   = y * grid_size.x() + x;
                        int idx = grid[i];
                        if (idx != -1) {
                            VectorDr p = gridPoints[idx];
                            if ((p - point).norm() < radius) { return false; }
                        }
                    }
                }
                return true;
            } else if constexpr (Dim == 3) {
                for (int z = std::max(0, grid_idx.z() - 2); z <= std::min(grid_size.z() - 1, grid_idx.z() + 2); ++z) {
                    for (int y = std::max(0, grid_idx.y() - 2); y <= std::min(grid_size.y() - 1, grid_idx.y() + 2);
                         ++y) {
                        for (int x = std::max(0, grid_idx.x() - 2); x <= std::min(grid_size.x() - 1, grid_idx.x() + 2);
                             ++x) {
                            int idx = grid[z * grid_size.y() * grid_size.x() + y * grid_size.x() + x];
                            if (idx != -1) {
                                VectorDr p = gridPoints[idx];
                                if ((p - point).norm() < radius) { return false; }
                            }
                        }
                    }
                }
                return true;
            }

            return false;
        }

        void addToGrid(const VectorDr & point) {
            VectorDi grid_idx = ((point - low_bound) / gridCellSize).cast<int>();
            if constexpr (Dim == 2) {
                grid[grid_idx.y() * grid_size.x() + grid_idx.x()] = gridPoints.size();
                gridPoints.push_back(point);
            } else if constexpr (Dim == 3) {
                grid[grid_idx.z() * grid_size.y() * grid_size.x() + grid_idx.y() * grid_size.x() + grid_idx.x()] =
                    gridPoints.size();
                gridPoints.push_back(point);
            }
        }
    };
} // namespace PhysX