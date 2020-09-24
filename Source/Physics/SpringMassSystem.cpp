#include "SpringMassSystem.h"

#include <iostream>
#include <cstdlib>

template <int Dim>
inline void SpringMassSystem<Dim>::Set_Particles(const std::vector<real> &_mass, const std::vector<VectorDr> &_positions, const std::vector<VectorDr> &_velocities)
{
	if (_mass.size() != _positions.size() || (!_velocities.empty() && _mass.size() != _velocities.size())) {
		std::cerr << "Error: [SpringMassSystem] set particles with mismatched sizes." << std::endl;
		std::exit(1);
	}
	num_particles = uint(_mass.size());
	num_dof = 3 * num_particles;
}

template <int Dim>
void SpringMassSystem<Dim>::Set_Springs(const std::vector<Spring> &_springs)
{
	springs = _springs;
}

template <int Dim>
void SpringMassSystem<Dim>::Set_Fixed_Particles(const std::vector<uint> &_fixed_particles)
{
	for (const uint fixed_particle : _fixed_particles) {
	}
}

template <int Dim>
void SpringMassSystem<Dim>::Write_Scene_Desc(std::ofstream &output) const
{
}

template <int Dim>
void SpringMassSystem<Dim>::Write_Frame(const std::string &frame_dir) const
{
}

template <int Dim>
void SpringMassSystem<Dim>::Save_Frame(const std::string &frame_dir) const
{
}

template <int Dim>
void SpringMassSystem<Dim>::Load_Frame(const std::string &frame_dir)
{
}

template <int Dim>
void SpringMassSystem<Dim>::Advance(const real dt)
{
}

template class SpringMassSystem<2>;
template class SpringMassSystem<3>;
