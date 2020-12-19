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

	void write(std::ofstream &fout) const;
	virtual void save(std::ofstream &fout) const;
	virtual void load(std::ifstream &fin);

	virtual void reinitialize();
	virtual void update(const real dt) = 0;
	virtual MatrixDr computeStressTensor(const int idx) = 0;
};

}
