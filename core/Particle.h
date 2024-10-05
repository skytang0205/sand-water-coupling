#pragma once

#include "Common.h"

namespace Pivot {
	struct Particle {
		Vector3d Position;
		Vector3d Velocity = Vector3d::Zero();
		std::array<Vector3d, 3> VelocityDrv = { Vector3d::Zero().eval(), Vector3d::Zero().eval(), Vector3d::Zero().eval() }; // Used for APIC
	};
}
