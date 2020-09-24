#pragma once

#include "Types.h"
#include "Simulation.h"

#include <string>

class Simulator
{
protected:

	Simulation *simulation;

	std::string output_dir;
	int first_frame;
	int last_frame;
	real frames_per_second;
	int steps_per_frame;

public:

	Simulator(
		Simulation *const _simulation,
		const std::string &_output_dir,
		const int _first_frame,
		const int _last_frame,
		const real _frames_per_second,
		const int _steps_per_frame
	) :
		simulation(_simulation),
		output_dir(_output_dir),
		first_frame(_first_frame),
		last_frame(_last_frame),
		frames_per_second(_frames_per_second),
		steps_per_frame(_steps_per_frame)
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
	real Get_Time_at_Frame(const int frame) const { return real(1) / frames_per_second * frame; }
};
