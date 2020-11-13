#pragma once

#include "Structures/GridBasedScalarField.h"
#include "Structures/GridBasedVectorField.h"

namespace PhysX {

template <int Dim>
class GridBasedAdvection
{
	DECLARE_DIM_TYPES(Dim)

public:

	GridBasedAdvection() = default;
	virtual ~GridBasedAdvection() = default;

	virtual void advect(FaceCenteredVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt) = 0;
};

template <int Dim>
class SemiLagrangianAdvection : public GridBasedAdvection<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	SemiLagrangianAdvection() = default;
	virtual ~SemiLagrangianAdvection() = default;

	void advect(FaceCenteredVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt);
};

}
