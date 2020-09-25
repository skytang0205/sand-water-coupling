#include "SpringMassSystem.h"
#include "File.h"

template <int Dim>
inline void SpringMassSystem<Dim>::Set_Particles(const std::vector<real> &_mass, const std::vector<VectorDr> &_positions, const std::vector<VectorDr> &_velocities)
{
	if (_mass.size() != _positions.size() || (!_velocities.empty() && _mass.size() != _velocities.size())) {
		std::cerr << "Error: [SpringMassSystem] set particles with mismatched sizes." << std::endl;
		std::exit(1);
	}

	num_particles = uint(_mass.size());
	num_dof = Dim * num_particles;
	
	mass.resize(num_dof);
	for (uint i = 0; i < num_particles; i++) {
		mass[Dim * i] = mass[Dim * i + 1] = _mass[i];
		if constexpr (Dim == 3) mass[Dim * i + 2] = _mass[i];
	}
	inv_mass = mass.cwiseInverse();
	inv_sqrt_mass = mass.cwiseSqrt().cwiseInverse();

	positions = _positions;
	velocities = _velocities;
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
		inv_mass[Dim * fixed_particle] = inv_mass[Dim * fixed_particle + 1] = 0;
		inv_sqrt_mass[Dim * fixed_particle] = inv_sqrt_mass[Dim * fixed_particle + 1] = 0;
		if constexpr (Dim == 3) inv_mass[Dim * fixed_particle + 2] = inv_sqrt_mass[Dim * fixed_particle + 2] = 0;
	}
}

template <int Dim>
void SpringMassSystem<Dim>::Write_Scene_Desc(std::ofstream &output) const
{
	output << "spring " << Dim << std::endl;
	output << "line_list indexed" << std::endl;
	output << "1.0 0.0 0.0 1.0" << std::endl;
	output << "0.01 0.01 0.01" << std::endl;
	output << "0.0" << std::endl;
}

template <int Dim>
void SpringMassSystem<Dim>::Write_Frame(const std::string &frame_dir) const
{
	// Write the springs.
	{
		std::string file_name = frame_dir + "/spring.mesh";
		std::ofstream output(file_name, std::ios::binary);
		File::Write_Binary(output, uint(num_particles));
		for (uint i = 0; i < num_particles; i++) {
			VectorDf pos = positions[i].cast<float>();
			Vector3f normal = Vector3f::Unit(2);
			File::Write_Binary(output, pos);
			File::Write_Binary(output, normal);
		}
		File::Write_Binary(output, uint(springs.size()) * 2);
		for (const auto &spring : springs) {
			File::Write_Binary(output, uint(spring.pid_1));
			File::Write_Binary(output, uint(spring.pid_2));
		}
	}
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
