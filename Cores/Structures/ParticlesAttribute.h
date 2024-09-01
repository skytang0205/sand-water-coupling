#pragma once

#include "Utilities/IO.h"
#include "Utilities/Types.h"

#include <algorithm>
#include <numeric>
#include <vector>

namespace PhysX {

    template<int Dim> class Particles;
    template<int Dim> class SmoothedParticles;
    template<int Dim> class BoundaryParticles;
    template<int Dim> class VirtualParticle;
    template<int Dim> class Shapes;

    template<int Dim, typename Type> class ParticlesAttribute {
        DECLARE_DIM_TYPES(Dim)

    public:
        friend class Particles<Dim>;
        friend class SmoothedParticles<Dim>;
        friend class BoundaryParticles<Dim>;
        friend class VirtualParticle<Dim>;
        friend class Shapes<Dim>;

    protected:
        std::vector<Type> _data;

    public:
        ParticlesAttribute()          = default;
        virtual ~ParticlesAttribute() = default;

        size_t size() const { return _data.size(); }
        bool   empty() const { return _data.empty(); }

        Type *       data() { return _data.data(); }
        const Type * data() const { return _data.data(); }

        Type &       operator[](const int idx) { return _data[idx]; }
        const Type & operator[](const int idx) const { return _data[idx]; }

        void setConstant(const Type & value) { std::fill(_data.begin(), _data.end(), value); }
        void setZero() { setConstant(Zero<Type>()); }

        template<typename AccType = Type> AccType sum() const {
            return std::accumulate(_data.begin(), _data.end(), Zero<AccType>());
        }

        Type min() const { return *std::min_element(_data.begin(), _data.end()); }
        Type max() const { return *std::max_element(_data.begin(), _data.end()); }

        Type absoluteMax() const {
            auto minmax = std::minmax_element(_data.begin(), _data.end());
            return std::max(std::abs(*minmax.first), std::abs(*minmax.second));
        }

        real normMax() const {
            if constexpr (HasSquaredNorm<Type>) {
                real squaredNormMax = 0;
                for (const auto & val : _data) { squaredNormMax = std::max(squaredNormMax, val.squaredNorm()); }
                return std::sqrt(squaredNormMax);
            } else return absoluteMax();
        }

        auto asVectorXr() {
            return Eigen::Map<VectorXr, Eigen::Aligned>(
                reinterpret_cast<real *>(_data.data()), _data.size() * (sizeof(Type) / sizeof(real)));
        }
        auto asVectorXr() const {
            return Eigen::Map<const VectorXr, Eigen::Aligned>(
                reinterpret_cast<const real *>(_data.data()), _data.size() * (sizeof(Type) / sizeof(real)));
        }

        void forEach(const std::function<void(const int)> & func) const {
            for (int i = 0; i < _data.size(); i++) func(i);
        }

        void parallelForEach(const std::function<void(const int)> & func) const {
#ifdef _OPENMP
#    pragma omp parallel for
#endif
            for (int i = 0; i < _data.size(); i++) func(i);
        }

        void load(std::istream & in) { IO::readArray(in, _data.data(), _data.size()); }
        void save(std::ostream & out) const { IO::writeArray(out, _data.data(), _data.size()); }
    };

    template<int Dim> using ParticlesScalarAttribute = ParticlesAttribute<Dim, real>;
    template<int Dim> using ParticlesVectorAttribute = ParticlesAttribute<Dim, Vector<Dim, real>>;

} // namespace PhysX
