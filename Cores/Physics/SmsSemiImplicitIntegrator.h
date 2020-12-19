#pragma once

#include "Physics/SpringMassSysIntegrator.h"
#include "Solvers/SparseSolver.h"

namespace PhysX {

template <int Dim>
class SmsSemiImplicitIntegrator : public SpringMassSysIntegrator<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using SpringMassSysIntegrator<Dim>::_particles;
	using SpringMassSysIntegrator<Dim>::_springs;

	ParticlesBasedVectorData<Dim> _forces;

	std::vector<Tripletr> _coeffBackwardEuler;
	SparseMatrixr _matBackwardEuler;
	VectorXr _rhsBackwardEuler;

	std::unique_ptr<SparseSolver> _solver;

public:

	SmsSemiImplicitIntegrator();

	SmsSemiImplicitIntegrator(const SmsSemiImplicitIntegrator &rhs) = delete;
	SmsSemiImplicitIntegrator &operator=(const SmsSemiImplicitIntegrator &rhs) = delete;
	virtual ~SmsSemiImplicitIntegrator() = default;

	virtual void reset(const Particles<Dim> *const particles, const std::vector<Spring> *const springs) override
	{
		SpringMassSysIntegrator<Dim>::reset(particles, springs);
		_forces.resize(particles);
		_matBackwardEuler.resize(particles->size() * Dim, particles->size() * Dim);
		_rhsBackwardEuler.resize(particles->size() * Dim);
	}

	virtual void integrate(
		const ParticlesVectorAttribute<Dim> &positions,
		ParticlesVectorAttribute<Dim> &velocities,
		const real dt,
		const ParticlesVectorAttribute<Dim> *const externalForces = nullptr,
		const std::unordered_set<int> *const constrainedDofs = nullptr) override;

protected:

	using SpringMassSysIntegrator<Dim>::accumulateForces;
};

}
