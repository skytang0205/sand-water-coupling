#pragma once

#include "Common.h"

namespace Pivot {
	struct Particle {
		Vector2d Position;
		Vector2d Velocity = Vector2d::Zero();
		std::array<Vector2d, 2> VelocityDrv = { Vector2d::Zero().eval(), Vector2d::Zero().eval() }; // Used for APIC
	};
}
