#pragma once

#include "Structures/Field.h"
#include "Structures/ParticlesBasedData.h"
#include "Structures/SmoothedParticles.h"

namespace PhysX {

template <int Dim>
class ParticlesBasedVectorField : public VectorField<Dim>, public ParticlesBasedVectorData<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using ParticlesVectorAttribute<Dim>::_data;
	using ParticlesBasedVectorData<Dim>::_particles;

public:

	ParticlesBasedVectorField(const SmoothedParticles<Dim> *const particles, const VectorDr &value = VectorDr::Zero()) :
		ParticlesBasedVectorData<Dim>(particles, value)
	{ }

	ParticlesBasedVectorField() = default;
	virtual ~ParticlesBasedVectorField() = default;

	void resize(const SmoothedParticles<Dim> *const particles, const VectorDr &value = VectorDr::Zero()) { ParticlesBasedVectorData<Dim>::resize(particles, value); }

	virtual VectorDr operator()(const VectorDr & pos) const override;
	virtual real divergence(const VectorDr & pos) const override;
};

}
