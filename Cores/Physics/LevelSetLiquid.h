#pragma once

#include "Geometries/GridBasedImplicitSurface.h"
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

	static constexpr int _kLsReinitMaxIters = 5;

	using EulerianFluid<Dim>::_kExtrapMaxIters;
	using EulerianFluid<Dim>::_grid;
	using EulerianFluid<Dim>::_velocity;
	using EulerianFluid<Dim>::_fluidFraction;
	using EulerianFluid<Dim>::_advector;
	using EulerianFluid<Dim>::_projector;

	GridBasedImplicitSurface<Dim> _levelSet;

	std::unique_ptr<LevelSetReinitializer<Dim>> _levelSetReinitializer;

	bool _enableGravity = true;
	
public:

	LevelSetLiquid(const StaggeredGrid<Dim> &grid);

	LevelSetLiquid(const LevelSetLiquid &rhs) = delete;
	LevelSetLiquid &operator=(const LevelSetLiquid &rhs) = delete;
	virtual ~LevelSetLiquid() = default;

	virtual void writeDescription(std::ofstream & fout) const override;
	virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string & frameDir) const override;
	virtual void loadFrame(const std::string & framdDir) override;

	virtual void initialize() override;

protected:

	using EulerianFluid<Dim>::enforceBoundaryConditions;

	virtual void advectFields(const real dt) override;
	virtual void applyBodyForces(const real dt) override;
	virtual void projectVelocity() override;

	virtual void extrapolateVelocity() override;
	virtual void reinitializeLevelSet();
};

}
