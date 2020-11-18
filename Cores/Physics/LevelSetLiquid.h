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

	using EulerianFluid<Dim>::_grid;
	using EulerianFluid<Dim>::_velocity;
	using EulerianFluid<Dim>::_advector;

	GridBasedImplicitSurface<Dim> _levelSet;

	std::unique_ptr<LevelSetReinitializer<Dim>> _levelSetReinitializer;
	
public:

	LevelSetLiquid(const StaggeredGrid<Dim> &grid);

	LevelSetLiquid(const LevelSetLiquid &rhs) = delete;
	LevelSetLiquid &operator=(const LevelSetLiquid &rhs) = delete;
	virtual ~LevelSetLiquid() = default;

	virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return real(1) / frameRate / stepRate; }

	virtual void writeDescription(std::ofstream & fout) const override;
	virtual void writeFrame(const std::string & frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string & frameDir) const override;
	virtual void loadFrame(const std::string & framdDir) override;

	virtual void initialize() override;
	virtual void advance(const real dt) override;
};

}
