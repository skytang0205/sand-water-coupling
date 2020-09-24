#pragma once

#include "Simulation.h"
#include "SurfaceMesh.h"

#include <vector>

template <int Dim>
class SpringMassSystem : public Simulation
{
private:

	DECLARE_DIM_TYPES(Dim)

public:

	enum struct Mode : uchar // name?
	{
		Euler,
		RK1 = Euler,
		Midpoint,
		RK2 = Midpoint,
		Heun,
		Modified_Euler = Heun,
		RK4,
		Backward_Euler,
		Trapezoidal_Rule
	};

	struct Spring
	{
		uint pid_1;
		uint pid_2;
		real rest_length;
		real stiffness;
		real damping;
	};

	friend class SpringMassSystemBuilder<Dim>;

protected:

	Mode mode = Mode::Backward_Euler;

	uint num_particles;
	uint num_dof;
	VectorXr mass;
	VectorXr inv_mass;
	VectorXr inv_sqrt_mass;
	std::vector<VectorDr> positions;
	std::vector<VectorDr> velocities;
	std::vector<Spring> springs;

	SurfaceMesh<Dim> surface_mesh; // std::unique_ptr?

	bool enable_gravity = false;

	SparseMatrixr A;

public:

	SpringMassSystem(const uint _num_particles) :
		num_particles(_num_particles),
		num_dof(Dim * num_particles),
		mass(num_dof),
		inv_mass(num_dof),
		inv_sqrt_mass(num_dof),
		positions(num_particles),
		velocities(num_particles),
		springs(0),
		surface_mesh(),
		A(num_dof, num_dof)
	{ }
};
