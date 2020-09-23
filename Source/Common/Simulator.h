#pragma once

#include "Types.h"
#include "Simulation.h"

#include <string>

struct SimulatorDesc
{
	std::string output_dir;
	int first_frame;
	int last_frame;
	real frames_per_second;
	int steps_per_frame;
};

class Simulator
{
protected:

	std::string output_dir;
	int first_frame;
	int last_frame;
	real frames_per_second;
	int steps_per_frame;

	Simulation *simulation;

public:

	Simulator(Simulation *const _simulation, const SimulatorDesc *const _simulator_desc) :
		simulation(_simulation),
		output_dir(_simulator_desc->output_dir),
		first_frame(_simulator_desc->first_frame),
		last_frame(_simulator_desc->last_frame),
		frames_per_second(_simulator_desc->frames_per_second),
		steps_per_frame(_simulator_desc->steps_per_frame)
	{ }

	Simulator(const Simulator &rhs) = delete;
	Simulator &operator=(const Simulator &rhs) = delete;
	virtual ~Simulator() = default;

	virtual void Simulate();

protected:

	void Create_Output_Directory() const;
	void Create_and_Write_Frame_Directory(const int frame) const;
	void Check_and_Load_Frame_Directory(const int frame);

	void Advance_Time_by_Steps(const real target_time);

	virtual real Get_Timestep() const { return real(1) / frames_per_second / steps_per_frame; }
	real Get_Time_at_Frame(const int frame) const { real(1) / frames_per_second * frame; }
};
