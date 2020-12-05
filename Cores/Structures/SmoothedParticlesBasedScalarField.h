#pragma once

#include "Structures/Field.h"
#include "Structures/ParticlesBasedData.h"
#include "Structures/SmoothedParticles.h"

namespace PhysX {

template <int Dim>
class SmoothedParticlesBasedScalarField : public ScalarField<Dim>, public ParticlesBasedScalarData<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

public:

};

}
