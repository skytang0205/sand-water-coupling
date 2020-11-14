#pragma once

#include "Structures/GridBasedScalarField.h"
#include "Structures/GridBasedVectorField.h"
#include "Structures/StaggeredGridBasedVectorField.h"

namespace PhysX {

template <int Dim>
class EulerianAdvector
{
	DECLARE_DIM_TYPES(Dim)

public:

	EulerianAdvector() = default;
	virtual ~EulerianAdvector() = default;

	virtual void advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt) = 0;
};

template <int Dim>
class SemiLagrangianAdvector : public EulerianAdvector<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	SemiLagrangianAdvector() = default;
	virtual ~SemiLagrangianAdvector() = default;

	void advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt);
};

}
