#pragma once

#include "Structures/Field.h"
#include "Structures/ParticlesBasedData.h"
#include "Structures/SmoothedParticles.h"

namespace PhysX {

    template<int Dim> class ParticlesBasedScalarField : public ScalarField<Dim>, public ParticlesBasedScalarData<Dim> {
        DECLARE_DIM_TYPES(Dim)

    protected:
        using ParticlesScalarAttribute<Dim>::_data;
        using ParticlesBasedScalarData<Dim>::_particles;

    public:
        ParticlesBasedScalarField(const SmoothedParticles<Dim> * const particles, const real value = 0):
            ParticlesBasedScalarData<Dim>(particles, value) {}

        ParticlesBasedScalarField()          = default;
        virtual ~ParticlesBasedScalarField() = default;

        void resize(const SmoothedParticles<Dim> * const particles, const real value = 0) {
            ParticlesBasedScalarData<Dim>::resize(particles, value);
        }

        virtual real     operator()(const VectorDr & pos) const override;
        virtual VectorDr gradient(const VectorDr & pos) const override;
        virtual real     laplacian(const VectorDr & pos) const override;

        VectorDr differenceGradientAtDataPoint(const int idx) const;
        VectorDr symmetricGradientAtDataPoint(const int idx) const;
        real     laplacianAtDataPoint(const int idx) const;
    };

}
