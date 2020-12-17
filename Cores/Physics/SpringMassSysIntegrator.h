#pragma once

#include "Solvers/SparseSolver.h"
#include "Structures/ParticlesBasedData.h"

#include <unordered_set>

namespace PhysX {

struct Spring
{
	int pid0;
	int pid1;
	real restLength;
	real stiffnessCoeff;
	real dampingCoeff;
};

template <int Dim>
class SpringMassSysIntegrator
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const Particles<Dim> *_particles = nullptr;
	const std::vector<Spring> *_springs = nullptr;

public:

	SpringMassSysIntegrator() = default;
	SpringMassSysIntegrator(const SpringMassSysIntegrator &rhs) = delete;
	SpringMassSysIntegrator &operator=(const SpringMassSysIntegrator &rhs) = delete;
	virtual ~SpringMassSysIntegrator() = default;

	virtual void reset(const Particles<Dim> *const particles, const std::vector<Spring> *const springs)
	{
		_particles = particles;
		_springs = springs;
	}

	virtual void integrate(
		const ParticlesVectorAttribute<Dim> &positions,
		ParticlesVectorAttribute<Dim> &velocities,
		const real dt,
		const ParticlesVectorAttribute<Dim> *const externalForces,
		const std::unordered_set<int> *const constrainedDofs) = 0;

protected:

	void accumulateAccelerations(
		const ParticlesVectorAttribute<Dim> &positions,
		const ParticlesVectorAttribute<Dim> &velocities,
		ParticlesVectorAttribute<Dim> &accelerations,
		const ParticlesVectorAttribute<Dim> *const externalForces = nullptr,
		const std::unordered_set<int> *const constrainedDofs = nullptr) const;
};

template <int Dim>
class SmsSymplecticEulerIntegrator : public SpringMassSysIntegrator<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using SpringMassSysIntegrator<Dim>::_particles;
	using SpringMassSysIntegrator<Dim>::_springs;

	ParticlesBasedVectorData<Dim> _accelerations;

public:

	SmsSymplecticEulerIntegrator() = default;
	SmsSymplecticEulerIntegrator(const SmsSymplecticEulerIntegrator &rhs) = delete;
	SmsSymplecticEulerIntegrator &operator=(const SmsSymplecticEulerIntegrator &rhs) = delete;
	virtual ~SmsSymplecticEulerIntegrator() = default;

	virtual void reset(const Particles<Dim> *const particles, const std::vector<Spring> *const springs) override
	{
		SpringMassSysIntegrator<Dim>::reset(particles, springs);
		_accelerations.resize(particles);
	}

	virtual void integrate(
		const ParticlesVectorAttribute<Dim> &positions,
		ParticlesVectorAttribute<Dim> &velocities,
		const real dt,
		const ParticlesVectorAttribute<Dim> *const externalForces = nullptr,
		const std::unordered_set<int> *const constrainedDofs = nullptr) override;

protected:

	using SpringMassSysIntegrator<Dim>::accumulateAccelerations;
};

template <int Dim>
class SmsBackwardEulerIntegrator : public SpringMassSysIntegrator<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using SpringMassSysIntegrator<Dim>::_particles;
	using SpringMassSysIntegrator<Dim>::_springs;

	ParticlesBasedVectorData<Dim> _accelerations;

	std::vector<Tripletr> _coeffBackwardEuler;
	SparseMatrixr _matBackwardEuler;
	VectorXr _rhsBackwardEuler;

	std::unique_ptr<SparseSolver> _solver;

public:

	SmsBackwardEulerIntegrator();

	SmsBackwardEulerIntegrator(const SmsBackwardEulerIntegrator &rhs) = delete;
	SmsBackwardEulerIntegrator &operator=(const SmsBackwardEulerIntegrator &rhs) = delete;
	virtual ~SmsBackwardEulerIntegrator() = default;

	virtual void reset(const Particles<Dim> *const particles, const std::vector<Spring> *const springs) override
	{
		SpringMassSysIntegrator<Dim>::reset(particles, springs);
		_accelerations.resize(particles);
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

	using SpringMassSysIntegrator<Dim>::accumulateAccelerations;
};

}
