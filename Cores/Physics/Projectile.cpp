#include "Projectile.h"

#include "Utilities/Constants.h"

#include <fstream>

namespace PhysX {

template<int Dim>
void Projectile<Dim>::writeDescription(std::ofstream &output) const
{
}

template<int Dim>
void Projectile<Dim>::writeFrame(const std::string &frameDir) const
{
	std::ofstream output(frameDir + "/ball.mesh", std::ios::binary);
	const VectorDf pos = _position.cast<float>();
	const float vel = float(_velocity.norm());
	output.write(reinterpret_cast<const char *>(&pos), sizeof(pos));
	output.write(reinterpret_cast<const char *>(&vel), sizeof(vel));
}

template<int Dim>
void Projectile<Dim>::saveFrame(const std::string &frameDir) const
{
	std::ofstream output(frameDir + "/state.sav", std::ios::binary);
	output.write(reinterpret_cast<const char *>(&_position), sizeof(_position));
	output.write(reinterpret_cast<const char *>(&_velocity), sizeof(_velocity));
}

template<int Dim>
void Projectile<Dim>::loadFrame(const std::string &frameDir)
{
	std::ifstream input(frameDir + "/state.sav", std::ios::binary);
	input.read(reinterpret_cast<char *>(&_position), sizeof(_position));
	input.read(reinterpret_cast<char *>(&_velocity), sizeof(_velocity));
}

template <int Dim>
void Projectile<Dim>::advance(const real dt)
{
	_position += _velocity * dt - VectorDr::Unit(1) * kGravity * dt * dt * real(0.5);
	_velocity.y() -= kGravity * dt;
	_time += dt;
}

template class Projectile<2>;
template class Projectile<3>;

}
