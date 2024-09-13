#pragma once

#include <algorithm>
#include <numeric>
#include <vector>
#include <random>
#include "Utilities/Types.h"
#include "Utilities/PoissonDiskSampler.h"
#include "Structures/ParticlesBasedData.h"
#include "Structures/ParticlesAttribute.h"
#include "Structures/StaggeredGrid.h"

namespace PhysX {

    template<int Dim> class Shapes{
        DECLARE_DIM_TYPES(Dim)
    public:
        ParticlesVectorAttribute<Dim> positions;
        ParticlesVectorAttribute<Dim> velocities;

        real _radius;

    public:
        Shapes(const real radius, const size_t cnt = 0, const VectorDr &pos = VectorDr::Zero(), const real mass = 1)
        {
            positions._data.resize(cnt, pos);
            velocities._data.resize(cnt, pos);
            _radius = radius;
        }
        virtual ~Shapes() = default;

        void generateBox(const VectorDr center = VectorDr::Zero(), const VectorDr & halfLengths = VectorDr::Ones(), bool if_Possion = false) {
            if (!if_Possion){
                if constexpr (Dim == 2) {
                    const real spacing = _radius * 2;
                    const int  xScale  = int(std::round(halfLengths[0] * 2 / spacing));
                    const int  yScale =
                        int(std::round((halfLengths[1] * 2 / spacing - 1) * 2 * real(std::numbers::inv_sqrt3) + 1));
                    const VectorDr origin =
                        center - VectorDr(xScale - 1, (yScale - real(1)) * real(std::numbers::sqrt3) / 2) * spacing / 2;
                    for (int j = 0; j < yScale; j++) {
                        const VectorDr yOrigin =
                            origin + VectorDr(j & 1 ? real(.5) : real(0), j * real(std::numbers::sqrt3) / 2) * spacing;
                        for (int i = 0; i < xScale - (j & 1); i++)
                            positions._data.push_back(yOrigin + VectorDr::Unit(0) * i * spacing);
                    }
                } else {
                    const real     spacing    = 2 * real(std::numbers::sqrt2) * _radius;
                    const VectorDi resolution = ((halfLengths - VectorDr::Ones() * _radius) * 2 / spacing)
                                                    .array()
                                                    .round()
                                                    .template cast<int>()
                                                    .matrix();
                    StaggeredGrid<Dim> grid(0, spacing, resolution, center);
                    grid.forEachFace(
                        [&](const int axis, const VectorDi & face) { positions._data.push_back(grid.faceCenter(axis, face)); });
                    grid.forEachNode([&](const VectorDi & node) { positions._data.push_back(grid.nodeCenter(node)); });
                }
            }
            else {
                // static std::random_device rd;
                // static std::mt19937 engine(rd());
                // static std::uniform_real_distribution dis(-1., 1.);

                // auto check_inside = [&](VectorDr pos) {
                //     for (const auto & p : positions._data)
                //         if ((p - pos).norm() < 2 * _radius)
                //             return true;
                //     return false;
                // };

                // if constexpr (Dim == 2) {
                //     int fail = 0;

                //     while (fail < 20000){
                //         double x = dis(engine), y = dis(engine);
                //         VectorDr rp = x * VectorDr::Unit(0) + y * VectorDr::Unit(1);
                //         VectorDr pos = center + rp * halfLengths.x();

                //         if (check_inside(pos))
                //             ++fail;
                //         else
                //             positions._data.push_back(pos), fail = 0;
                //     }
                // }

                PoissonDiskSampler<Dim> sampler(halfLengths, center, 2*_radius);
                std::vector<VectorDr> newPoints = sampler.generatePoints();
                positions._data.insert(positions._data.end(), newPoints.begin(), newPoints.end());
            }

            velocities._data.resize(positions.size());
        }

        void generateRotate(const real omega, const VectorDr &center = VectorDr::Zero()){
            velocities._data.resize(positions.size());
            velocities.parallelForEach([&](const int i) {
                velocities[i] = omega * (positions[i].y() - center.y()) * VectorDr::Unit(0)
                    - omega * (positions[i].x() - center.x()) * VectorDr::Unit(1);
            });
        }
    };
}