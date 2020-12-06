#pragma once

#include "Structures/Particles.h"

namespace PhysX {

template <int Dim, typename Type>
class ParticlesBasedData : public ParticlesAttribute<Dim, Type>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using ParticlesAttribute<Dim, Type>::_data;

	const Particles<Dim> *_particles = nullptr;

public:

	ParticlesBasedData(const Particles<Dim> *const particles, const Type &value = Zero<Type>()) { resize(particles, value); }

	ParticlesBasedData() = default;
	virtual ~ParticlesBasedData() = default;

	void resize(const Particles<Dim> *const particles, const Type &value = Zero<Type>())
	{
		_particles = particles;
		_data.resize(_particles->size(), value);
	}
};

template <int Dim> using ParticlesBasedScalarData = ParticlesBasedData<Dim, real>;
template <int Dim> using ParticlesBasedVectorData = ParticlesBasedData<Dim, Vector<Dim, real>>;

}
