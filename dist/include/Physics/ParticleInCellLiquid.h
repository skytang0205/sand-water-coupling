#pragma once

#include "Physics/LevelSetLiquid.h"
#include "Structures/ParticlesBasedData.h"

namespace PhysX {

template <int Dim>
class ParticleInCellLiquid : public LevelSetLiquid<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	friend class ParticleInCellLiquidBuilder;

protected:

	using EulerianFluid<Dim>::_kExtrapMaxSteps;
	using EulerianFluid<Dim>::_grid;
	using EulerianFluid<Dim>::_velocity;
	using EulerianFluid<Dim>::_advector;
	using EulerianFluid<Dim>::_boundaryHelper;
	using LevelSetLiquid<Dim>::_kLsReinitMaxSteps;
	using LevelSetLiquid<Dim>::_levelSet;
	using LevelSetLiquid<Dim>::_levelSetReinitializer;

	const int _particlesCntPerSubCell;
	Particles<Dim> _particles;
	ParticlesBasedVectorData<Dim> _particleVelocities;

public:

	ParticleInCellLiquid(const StaggeredGrid<Dim> &grid, const int particlesCntPerSubcell);

	ParticleInCellLiquid(const ParticleInCellLiquid &rhs) = delete;
	ParticleInCellLiquid &operator=(const ParticleInCellLiquid &rhs) = delete;
	virtual ~ParticleInCellLiquid() = default;

	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

	virtual void initialize() override;
	virtual void advance(const real dt) override;

protected:

	using EulerianFluid<Dim>::updateColliders;
	using LevelSetLiquid<Dim>::applyBodyForces;
	using LevelSetLiquid<Dim>::projectVelocity;

	virtual void advectFields(const real dt) override;
	virtual void applyParticleForces(const real dt) { }
	virtual void transferFromGridToParticles();
	virtual void transferFromParticlesToGrid();

	virtual void transferFromParticlesToGrid(StaggeredGridBasedScalarData<Dim> &weightSum);
	virtual void maintainGridBasedData(StaggeredGridBasedScalarData<Dim> &weightSum);
	virtual void reinitializeLevelSet() override;
	virtual void reinitializeParticles();
	virtual void reinitializeParticlesBasedData();
};

}
