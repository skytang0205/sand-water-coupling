#pragma once

#include "Utilities/Types.h"

namespace PhysX {

struct Spring
{
	int pid0;
	int pid1;
	real restLength;
	real stiffnessCoeff;
	real dampingCoeff;
};

}
