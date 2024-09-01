#pragma once

#include "Structures/GridBasedScalarField.h"
#include "Structures/GridBasedVectorField.h"
#include "Structures/ParticlesAttribute.h"
#include "Structures/StaggeredGridBasedVectorField.h"

namespace PhysX {

template <int Dim>
class EulerianAdvector
{
	DECLARE_DIM_TYPES(Dim)

public:

	EulerianAdvector() = default;
	EulerianAdvector(const EulerianAdvector &rhs) = delete;
	EulerianAdvector &operator=(const EulerianAdvector &rhs) = delete;
	virtual ~EulerianAdvector() = default;

	virtual void advect(GridBasedScalarField<Dim> &field, const VectorField<Dim> &flow, const real dt) = 0;
	virtual void advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt) = 0;
	virtual void advect(ParticlesVectorAttribute<Dim> &positions, const VectorField<Dim> &flow, const real dt) = 0;
};

template <int Dim, int RungeKuttaOrder = 2>
class SemiLagrangianAdvector : public EulerianAdvector<Dim>
{
	DECLARE_DIM_TYPES(Dim)

	static_assert(1 <= RungeKuttaOrder && RungeKuttaOrder <= 4, "Runge-Kutta order must be 1, 2, 3 or 4.");

public:

	SemiLagrangianAdvector() = default;
	SemiLagrangianAdvector(const SemiLagrangianAdvector &rhs) = delete;
	SemiLagrangianAdvector &operator=(const SemiLagrangianAdvector &rhs) = delete;
	virtual ~SemiLagrangianAdvector() = default;

	virtual void advect(GridBasedScalarField<Dim> &field, const VectorField<Dim> &flow, const real dt) override;
	virtual void advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt) override;
	virtual void advect(ParticlesVectorAttribute<Dim> &positions, const VectorField<Dim> &flow, const real dt) override;

protected:

	void advect(const GridBasedScalarField<Dim> &field, GridBasedScalarField<Dim> &newField, const VectorField<Dim> &flow, const real dt) const;
	void advect(const StaggeredGridBasedVectorField<Dim> &field, StaggeredGridBasedVectorField<Dim> &newField, const VectorField<Dim> &flow, const real dt) const;

	VectorDr trace(const VectorDr &startPos, const VectorField<Dim> &flow, const real dt) const;
};

template <int Dim, int RungeKuttaOrder = 2>
class MacCormackAdvector : public SemiLagrangianAdvector<Dim, RungeKuttaOrder>
{
	DECLARE_DIM_TYPES(Dim)

	static_assert(1 <= RungeKuttaOrder && RungeKuttaOrder <= 4, "Runge-Kutta order must be 1, 2, 3 or 4.");

public:

	MacCormackAdvector() = default;
	MacCormackAdvector(const MacCormackAdvector & rhs) = delete;
	MacCormackAdvector &operator=(const MacCormackAdvector & rhs) = delete;
	virtual ~MacCormackAdvector() = default;

	virtual void advect(GridBasedScalarField<Dim> &field, const VectorField<Dim> &flow, const real dt) override;
	virtual void advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt) override;
	virtual void advect(ParticlesVectorAttribute<Dim> &positions, const VectorField<Dim> &flow, const real dt) override;

protected:

	using SemiLagrangianAdvector<Dim, RungeKuttaOrder>::trace;
};

}
