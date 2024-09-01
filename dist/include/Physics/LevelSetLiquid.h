#pragma once

#include "Geometries/LevelSetReinitializer.h"
#include "Physics/EulerianFluid.h"

namespace PhysX {

template <int Dim>
class LevelSetLiquid : public EulerianFluid<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	friend class LevelSetLiquidBuilder;

protected:

	static constexpr int _kLsReinitMaxSteps = 5;

	using EulerianFluid<Dim>::_kExtrapMaxSteps;
	using EulerianFluid<Dim>::_grid;
	using EulerianFluid<Dim>::_velocity;
	using EulerianFluid<Dim>::_colliders;
	using EulerianFluid<Dim>::_advector;
	using EulerianFluid<Dim>::_boundaryHelper;
	using EulerianFluid<Dim>::_projector;

	LevelSet<Dim> _levelSet;

	std::unique_ptr<LevelSetReinitializer<Dim>> _levelSetReinitializer;

	bool _enableGravity = true;
	bool _enableSurfaceTension = false;

	real _density = real(1e3);
	real _surfaceTensionCoefficient = real(7.28e-2);

public:

	LevelSetLiquid(const StaggeredGrid<Dim> &grid);

	LevelSetLiquid(const LevelSetLiquid &rhs) = delete;
	LevelSetLiquid &operator=(const LevelSetLiquid &rhs) = delete;
	virtual ~LevelSetLiquid() = default;

	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string & frameDir) const override;
	virtual void loadFrame(const std::string & frameDir) override;

	virtual void initialize() override;

protected:

	virtual void advectFields(const real dt) override;
	virtual void applyBodyForces(const real dt) override;
	virtual void projectVelocity(const real dt = 0) override;

	virtual void reinitializeLevelSet();
};

}
