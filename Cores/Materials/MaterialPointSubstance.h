#pragma once

#include "Structures/ParticlesBasedData.h"
#include "Utilities/IO.h"

namespace PhysX {

template <int Dim>
class MaterialPointSubstance
{
	DECLARE_DIM_TYPES(Dim)

public:

	Particles<Dim> particles;
	ParticlesBasedVectorData<Dim> velocities;
	ParticlesBasedData<Dim, MatrixDr> velocityDerivatives;

protected:

	const std::string _name;
	const Vector4f _color;

	const real _density;

public:

	MaterialPointSubstance(const std::string &name, const Vector4f &color, const real density) : _name(name), _color(color), _density(density) { }

	virtual ~MaterialPointSubstance() = default;

	const std::string &name() const { return _name; }
	const Vector4f &color() const { return _color; }
	real density() const { return _density; }

	void write(std::ofstream &fout) const
	{
		IO::writeValue(fout, uint(particles.size()));
		particles.forEach([&](const int i) {
			IO::writeValue(fout, particles.positions[i].template cast<float>().eval());
		});
		if constexpr (Dim == 3) {
			particles.forEach([&](const int i) {
				IO::writeValue(fout, VectorDf::Unit(2).eval());
			});
		}
	}

	virtual void save(std::ofstream &fout) const
	{
		IO::writeValue(fout, uint(particles.size()));
		particles.positions.save(fout);
	}

	virtual void load(std::ifstream &fin)
	{
		uint particlesCnt;
		IO::readValue(fin, particlesCnt);
		particles.resize(particlesCnt);
		particles.positions.load(fin);
		reinitialize();
	}

	virtual void reinitialize()
	{
		velocities.resize(&particles);
		velocities.setZero();

		velocityDerivatives.resize(&particles);
		velocityDerivatives.setZero();
	}

	virtual void update(const int idx, const real dt) { particles.positions[idx] += velocities[idx] * dt; }
	virtual MatrixDr computeStressTensor(const int idx) const = 0;
};

}
