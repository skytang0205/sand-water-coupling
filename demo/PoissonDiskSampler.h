#pragma once

#include "Surface.h"
#include "Common.h"
#include "ImplicitSurface.h"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <vector>

namespace Pivot {
    class PoissonDiskSampler {
    public:
        PoissonDiskSampler(double radius, Surface const &surface):
            radius(radius){
            auto const [minCorner, maxCorner] = surface.GetCornersOfAABB();
            low_bound = minCorner;
            up_bound = maxCorner;
            gridCellSize = radius / std::numbers::sqrt3;
            grid_size    = ((up_bound - low_bound) / gridCellSize).cast<int>();
            int size     = 1;
            for (int i = 0; i < 3; i++) size *= grid_size(i);
            grid.resize(size, -1);
        }

        Vector3d generateNewPoint(Vector3d point) {
            static std::random_device               rd;
            static std::mt19937                     gen(rd());
            static std::uniform_real_distribution<> dis(0.0, 1.0);

            double            angle1   = dis(gen) * 2 * std::numbers::pi;
            double            angle2   = std::acos(1 - 2 * dis(gen));
            double            dist     = radius + dis(gen) * radius;
            Vector3d newPoint = { std::cos(angle1) * std::sin(angle2),
                                            std::sin(angle1) * std::sin(angle2),
                                            std::cos(angle2) };

            newPoint = newPoint * dist + point;
            return newPoint;
            
        }

        std::vector<Vector3d> generatePoints(Surface const &surface) {
            std::vector<Vector3d> samples;
            std::vector<Vector3d> processList;
            Vector3d              initialPoint = (up_bound + low_bound) * 0.5;
            processList.push_back(initialPoint);
            addToGrid(initialPoint);

            while (! processList.empty()) {
                int      idx   = rand() % processList.size();
                Vector3d point = processList[idx];
                bool     found = false;

                for (int i = 0; i < 3000; ++i) {
                    Vector3d newPoint = generateNewPoint(point);
                    if (isInBounds(newPoint) && surface.Surrounds(newPoint) && isFarEnough(newPoint)) {
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
                if (gridPoints.size() > grid.size() * 10) break;
            }

            for (const auto & idx : grid) {
                if (idx != -1) { samples.push_back(gridPoints[idx]); }
            }
            return samples;
        }

        bool isInBounds(const Vector3d & point) const {
            Vector3i grid_idx = ((point - low_bound) / gridCellSize).cast<int>();
            return grid_idx.x() >= 0 && grid_idx.x() < grid_size.x() && grid_idx.y() >= 0
                    && grid_idx.y() < grid_size.y() && grid_idx.z() >= 0 && grid_idx.z() < grid_size.z();
        }

    private:
        double                radius, gridCellSize;
        Vector3d              low_bound, up_bound;
        Vector3i              grid_size;
        std::vector<int>      grid;
        std::vector<Vector3d> gridPoints;

        bool isFarEnough(const Vector3d & point) const {
            Vector3i grid_idx = ((point - low_bound) / gridCellSize).cast<int>();

            // for (const auto & p : gridPoints)
            //     if ((p - point).norm() < radius)
            //         return false;
            // return true;

            for (int z = std::max(0, grid_idx.z() - 2); z <= std::min(grid_size.z() - 1, grid_idx.z() + 2); ++z) {
                for (int y = std::max(0, grid_idx.y() - 2); y <= std::min(grid_size.y() - 1, grid_idx.y() + 2);
                        ++y) {
                    for (int x = std::max(0, grid_idx.x() - 2); x <= std::min(grid_size.x() - 1, grid_idx.x() + 2);
                            ++x) {
                        int idx = grid[z * grid_size.y() * grid_size.x() + y * grid_size.x() + x];
                        if (idx != -1) {
                            Vector3d p = gridPoints[idx];
                            if ((p - point).norm() < radius) { return false; }
                        }
                    }
                }
            }
            return true;
        }

        void addToGrid(const Vector3d & point) {
            Vector3i grid_idx = ((point - low_bound) / gridCellSize).cast<int>();
            grid[grid_idx.z() * grid_size.y() * grid_size.x() + grid_idx.y() * grid_size.x() + grid_idx.x()] =
                gridPoints.size();
            gridPoints.push_back(point);
            
        }
    };
} // namespace PhysX