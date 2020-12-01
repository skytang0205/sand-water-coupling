#pragma once

#include "Geometries/Collider.h"
#include "Physics/EulerianAdvector.h"
#include "Physics/EulerianBoundaryHelper.h"
#include "Physics/EulerianProjector.h"
#include "Physics/Simulation.h"
#include "Structures/StaggeredGridBasedData.h"
#include "Structures/StaggeredGridBasedVectorField.h"

#include <functional>
#include <memory>

namespace PhysX {

template <int Dim>
class EulerianFluid : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

public:

	friend class EulerianFluidBuilder;

protected:

	static constexpr int _kExtrapMaxSteps = 3;

	const StaggeredGrid<Dim> _grid;

	StaggeredGridBasedVectorField<Dim> _velocity;

	std::vector<std::unique_ptr<Collider<Dim>>> _colliders;
	std::function<real(const int, const VectorDi &)> _domainBoundaryVelocity;

	std::unique_ptr<EulerianBoundaryHelper<Dim>> _boundaryHelper;
	std::unique_ptr<EulerianAdvector<Dim>> _advector;
	std::unique_ptr<EulerianProjector<Dim>> _projector;

public:

	EulerianFluid(const StaggeredGrid<Dim> &grid);

	EulerianFluid(const EulerianFluid &rhs) = delete;
	EulerianFluid &operator=(const EulerianFluid &rhs) = delete;
	virtual ~EulerianFluid() = default;

	virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return std::min(real(1) / frameRate, stepRate * _grid.spacing() / _velocity.absoluteMax()); }

	virtual int dimension() const override { return Dim; }
	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

	virtual void initialize() override;
	virtual void advance(const real dt) override;

protected:

	virtual void advectFields(const real dt);
	virtual void updateColliders(const real dt);
	virtual void applyBodyForces(const real dt) { }
	virtual void projectVelocity(const real dt = 0);

	virtual void updateBoundary();
};

}
