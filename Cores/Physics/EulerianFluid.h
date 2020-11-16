#pragma once

#include "Geometries/Collider.h"
#include "Geometries/GridBasedImplicitSurface.h"
#include "Physics/EulerianAdvector.h"
#include "Physics/EulerianProjector.h"
#include "Physics/Simulation.h"
#include "Structures/StaggeredGridBasedData.h"
#include "Structures/StaggeredGridBasedVectorField.h"

#include <memory>

namespace PhysX {

template <int Dim>
class EulerianFluid : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const StaggeredGrid<Dim> _grid;

	StaggeredGridBasedVectorField<Dim> _velocity;
	StaggeredGridBasedData<Dim> _fluidFraction;

	std::vector<std::unique_ptr<Collider<Dim>>> _colliders;

	std::unique_ptr<EulerianAdvector<Dim>> _advector;
	std::unique_ptr<EulerianProjector<Dim>> _projector;

public:

	EulerianFluid(const StaggeredGrid<Dim> &grid);

	EulerianFluid(const EulerianFluid &rhs) = delete;
	EulerianFluid &operator=(const EulerianFluid &rhs) = delete;
	virtual ~EulerianFluid() = default;

	virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return real(1) / frameRate / stepRate; }

	virtual void writeDescription(std::ofstream &fout) const;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const;
	virtual void saveFrame(const std::string &frameDir) const;
	virtual void loadFrame(const std::string &framdDir);

	virtual void initialize();
	virtual void advance(const real dt);

protected:

	void advectFields(const real dt);
	void updateColliders(const real dt);
	void applyBodyForces(const real dt);
	void projectVelocity();

	void updateFluidFraction();
	void extrapolateVelocity(const int maxIterations = 3);
	void enforceBoundaryConditions();
};

}
