#pragma once

#include "Structures/Particles.h"

namespace PhysX {

template <int Dim, typename Type>
class ParticlesBasedData
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const Particles<Dim> *_particles = nullptr;
	std::vector<Type> _data;

public:

	ParticlesBasedData(const Particles<Dim> *const particles, const Type &value = Zero<Type>()) { resize(particles, value); }

	ParticlesBasedData() = default;
	virtual ~ParticlesBasedData() = default;

	void resize(const Particles<Dim> *const particles, const Type &value = Zero<Type>())
	{
		_particles = particles;
		_data.resize(_particles->size(), value);
	}

	size_t size() const { return _particles->size(); }
	bool empty() const { return _particles->empty(); }

	Type *data() { return _data.data(); }
	const Type *data() const { return _data.data(); }

	Type &operator[](const size_t index) { return _data[index]; }
	const Type &operator[](const size_t index) const { return _data[index]; }

	void setConstant(const Type &value) { std::fill(_data.begin(), _data.end(), value); }
	void setZero() { setConstant(Zero<Type>()); }

	void forEach(const std::function<void(const int)> &func) const { _particles->forEach(func); }
	void parallelForEach(const std::function<void(const int)> &func) const { _particles->parallelForEach(func); }

	void load(std::istream &in) { IO::readArray(in, _data.data(), _data.size()); }
	void save(std::ostream &out) const { IO::writeArray(out, _data.data(), _data.size()); }
};

template <int Dim> using ParticlesBasedScalarData = ParticlesBasedData<Dim, real>;
template <int Dim> using ParticlesBasedVectorData = ParticlesBasedData<Dim, Vector<Dim, real>>;

}
