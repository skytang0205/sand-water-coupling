#pragma once

#include "Physics/Simulation.h"
#include "Solvers/SparseSolver.h"
#include "Structures/ParticlesAttribute.h"

namespace PhysX {

template <int Dim>
class SpringMassSystem : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

public:

	friend class SpringMassSystemBuilder;

	struct Spring
	{
		int pid0;
		int pid1;
		real restLength;
		real stiffnessCoeff;
		real dampingCoeff;
	};

protected:

	ParticlesVectorAttribute<Dim> _masses;
	ParticlesVectorAttribute<Dim> _invMasses;
	ParticlesVectorAttribute<Dim> _invSqrtMasses;
	ParticlesVectorAttribute<Dim> _positions;
	ParticlesVectorAttribute<Dim> _velocities;
	ParticlesVectorAttribute<Dim> _accelerations;
	std::vector<Spring> _springs;

	std::vector<Tripletr> _coeffBackwardEuler;
	SparseMatrixr _matBackwardEuler;
	VectorXr _rhsBackwardEuler;

	std::unique_ptr<SparseSolver> _solver;

	bool _enableGravity = true;

public:

	SpringMassSystem();

	SpringMassSystem(const SpringMassSystem &rhs) = delete;
	SpringMassSystem &operator=(const SpringMassSystem &rhs) = delete;
	virtual ~SpringMassSystem() = default;

	virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return real(1) / frameRate / stepRate; }

	virtual int dimension() const override { return Dim; }
	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

	virtual void initialize() override;
	virtual void advance(const real dt) override;

protected:

	void calculateAccelarations();
	void reinitializeAttributes();
	void buildAndSolveLinearSystem(const real dt);

	virtual void applyElasticForces();
	virtual void applyExternalForces();
};

}
