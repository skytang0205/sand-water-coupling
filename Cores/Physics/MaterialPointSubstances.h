#pragma once

#include "Geometries/Collider.h"
#include "Physics/Simulation.h"
#include "Structures/GridBasedData.h"
#include "Structures/ParticlesBasedData.h"
#include "Structures/StaggeredGrid.h"

#include <string>

namespace PhysX {

template <int Dim>
class MaterialPointSubstances : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

public:

	friend class MatPointSubstancesBuilder;

	class Substance
	{
	public:

		Particles<Dim> particles;
		ParticlesBasedVectorData<Dim> velocities;
		ParticlesBasedData<Dim, MatrixDr> velocityDerivatives;
		ParticlesBasedData<Dim, MatrixDr> deformationGradients;
		ParticlesBasedScalarData<Dim> plasticJacobians;

	protected:

		const std::string _name;
		const Vector4f _color;

		const real _density;
		const real _lameLambda;
		const real _lameMu;

		const real _plasticLowerBound;
		const real _plasticUpperBound;
		const real _hardeningCoeff;

	public:

		Substance(const std::string &name, const Vector4f &color, const real density, const real lambda, const real mu, const real thetaC = -1, const real thetaS = -1, const real hardeningCoeff = 0) :
			_name(name),
			_color(color),
			_density(density),
			_lameLambda(lambda),
			_lameMu(mu),
			_plasticLowerBound(1 - thetaC),
			_plasticUpperBound(1 + thetaS),
			_hardeningCoeff(hardeningCoeff)
		{ }

		virtual ~Substance() = default;

		const std::string &name() const { return _name; }
		const Vector4f &color() const { return _color; }
		real density() const { return _density; }
		bool plastic() const { return _plasticLowerBound <= _plasticUpperBound; }

		MatrixDr computeStressTensor(const int idx);
		void reinitialize();
	};

protected:

	std::vector<Substance> _substances;

	const StaggeredGrid<Dim> _grid;

	GridBasedVectorData<Dim> _velocity;
	GridBasedScalarData<Dim> _mass;

	const StaticCollider<Dim> _domainBoundary;
	std::vector<std::unique_ptr<Collider<Dim>>> _colliders;

	bool _enableGravity = true;

public:

	MaterialPointSubstances(const StaggeredGrid<Dim> &grid);

	MaterialPointSubstances(const MaterialPointSubstances &rhs) = delete;
	MaterialPointSubstances &operator=(const MaterialPointSubstances &rhs) = delete;
	virtual ~MaterialPointSubstances() = default;

	virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return real(1) / frameRate / stepRate; }

	virtual int dimension() const override { return Dim; }
	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

	virtual void initialize() override;
	virtual void advance(const real dt) override;

protected:

	virtual void applyLagrangianForces(const real dt) { }
	virtual void applyEulerianForces(const real dt);

	virtual void transferFromParticlesToGrid(const real dt);
	virtual void transferFromGridToParticles(const real dt);

	void sampleParticlesInsideSurface(Substance &substance, const Surface<Dim> &surface, const int particlesCntPerSubcell);
};

}
