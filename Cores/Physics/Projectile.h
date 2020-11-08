#pragma once

#include "Physics/Simulation.h"

namespace PhysX {

template <int Dim>
class Projectile : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

protected:

	VectorDr _position;
	VectorDr _velocity;

public:

	Projectile(const VectorDr &position, const VectorDr &velocity) :
		_position(position),
		_velocity(velocity)
	{ }

	Projectile() = delete;
	Projectile(const Projectile &rhs) = delete;
	Projectile &operator=(const Projectile &rhs) = delete;
	virtual ~Projectile() = default;

	virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return real(1) / frameRate / stepRate; }

	virtual void writeDescription(std::ofstream &output) const;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw = false) const;
	virtual void saveFrame(const std::string &frameDir) const;
	virtual void loadFrame(const std::string &framdDir);

	virtual void advance(const real dt);
};

}
