#pragma once

#include "Physics/GirdBasedAdvection.h"
#include "Physics/Simulation.h"
#include "Structures/GridBasedVectorField.h"

#include <memory>

namespace PhysX {

template <int Dim>
class GridBasedFluid : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const StaggeredGrid<Dim> _grid;

	GridBasedStaggeredVectorField<Dim> _velocity;

	std::unique_ptr<GridBasedAdvection<Dim>> _advection;

public:

	GridBasedFluid(const StaggeredGrid<Dim> &grid);

	GridBasedFluid(const GridBasedFluid &rhs) = delete;
	GridBasedFluid &operator=(const GridBasedFluid &rhs) = delete;
	virtual ~GridBasedFluid() = default;

	virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return real(1) / frameRate / stepRate; }

	virtual void writeDescription(std::ofstream &fout) const;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const;
	virtual void saveFrame(const std::string &frameDir) const;
	virtual void loadFrame(const std::string &framdDir);

	virtual void advance(const real dt);
};

}
