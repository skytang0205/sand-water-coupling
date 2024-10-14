#pragma once

#include "Common.h"

namespace Pivot {
	enum class ParticleType : std::uint32_t {
		PIC  = 0,
		DEM  = 1,
	};

	struct Particle {
		ParticleType Type = ParticleType::PIC;
		Vector2d Position;
		Vector2d Velocity = Vector2d::Zero();
		Vector2d AccVelocity = Vector2d::Zero();
		Vector2d CouplingForce = Vector2d::Zero();
		std::array<Vector2d, 2> VelocityDrv = { Vector2d::Zero().eval(), Vector2d::Zero().eval() }; // Used for APIC
	};
}
