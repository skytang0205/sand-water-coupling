#pragma once

#include "Geometries/Collider.h"
#include "Physics/Simulation.h"
#include "Structures/GridBasedData.h"
#include "Structures/ParticlesBasedData.h"
#include "Structures/StaggeredGrid.h"

namespace PhysX {

template <int Dim>
class MaterialPointSubstances : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

public:

	struct Substance
	{
		Particles<Dim> particles;
		ParticlesBasedVectorData<Dim> velocities;
		ParticlesBasedData<Dim, MatrixDr> velocityDerivatives;
		ParticlesBasedData<Dim, MatrixDr> deformationGradients;
		ParticlesBasedScalarData<Dim> volumeRatios;

		Vector4f color;
	};

protected:

	std::vector<Substance> _substances;
	GridBasedVectorData<Dim> _velocity;
	GridBasedScalarData<Dim> _mass;

	std::vector<std::unique_ptr<Collider<Dim>>> _colliders;

public:

	MaterialPointSubstances(const StaggeredGrid<Dim> &grid);

	MaterialPointSubstances(const MaterialPointSubstances &rhs) = delete;
	MaterialPointSubstances &operator=(const MaterialPointSubstances &rhs) = delete;
	virtual ~MaterialPointSubstances() = default;

	virtual int dimension() const override { return Dim; }
	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

	virtual void initialize() override;
	virtual void advance(const real dt) override;

protected:

	virtual void advect(const real dt);

	virtual void transferFromGridToParticles();
	virtual void transferFromParticlesToGrid();
};

}
