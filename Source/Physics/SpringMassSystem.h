#pragma once

#include "Simulation.h"
#include "SurfaceMesh.h"

#include <vector>

template <int Dim> class SpringMassSystemBuilder;

template <int Dim>
class SpringMassSystem : public Simulation
{
private:

	DECLARE_DIM_TYPES(Dim)

public:

	enum struct IntegratorType : uchar
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

	IntegratorType integrator_type = IntegratorType::Backward_Euler;

	uint num_particles = 0;
	uint num_dof = 0;
	VectorXr mass;
	VectorXr inv_mass;
	VectorXr inv_sqrt_mass;
	std::vector<VectorDr> positions;
	std::vector<VectorDr> velocities;
	std::vector<Spring> springs;

	bool enable_gravity = false;

	SparseMatrixr A;

public:

	SpringMassSystem() = default;
	SpringMassSystem(const SpringMassSystem &rhs) = delete;
	SpringMassSystem &operator=(const SpringMassSystem &rhs) = delete;
	virtual ~SpringMassSystem() = default;

	void Set_Particles(const std::vector<real> &_mass, const std::vector<VectorDr> &_positions, const std::vector<VectorDr> &_velocities = {});
	void Set_Springs(const std::vector<Spring> &_springs);
	void Set_Fixed_Particles(const std::vector<uint> &_fixed_particles);

	virtual void Write_Scene_Desc(std::ofstream &output) const override;
	virtual void Write_Frame(const std::string &frame_dir) const override;
	virtual void Save_Frame(const std::string &frame_dir) const override;
	virtual void Load_Frame(const std::string &frame_dir) override;

	virtual void Advance(const real dt) override;
};
