#pragma once

#include "Geometries/Collider.h"
#include "Physics/Simulation.h"
#include "Physics/SpringMassSysIntegrator.h"
#include "Structures/ParticlesBasedData.h"

#include <unordered_set>

namespace PhysX {

template <int Dim>
class SpringMassSystem : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

public:

	friend class SpringMassSystemBuilder;

protected:

	Particles<Dim> _particles;
	ParticlesBasedVectorData<Dim> _velocities;
	ParticlesBasedVectorData<Dim> _externalForces;

	std::vector<Spring> _springs;
	std::unordered_set<int> _constrainedDofs;

	std::vector<std::unique_ptr<Collider<Dim>>> _colliders;

	std::unique_ptr<SpringMassSysIntegrator<Dim>> _integrator;

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

	void reinitializeParticlesBasedData();
	void moveParticles(const real dt);

	void calculateExternalForces();
	void resolveVelocities(const real dt);
};

}
