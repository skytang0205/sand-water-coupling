#pragma once

#include "Types.h"

#include <string>

// the base class of various simulations
class Simulation
{
protected:

	real time; // the current time

public:

	Simulation(const real _time = 0) : time(_time) { }

	Simulation(const Simulation &rhs) = delete;
	Simulation &operator=(const Simulation &rhs) = delete;
	virtual ~Simulation() = default;

	void Set_Time(const real new_time = 0) { time = new_time; }
	real Get_Time() const { return time; }

	virtual void Write_Scene_Desc(std::ofstream &output) const = 0;
	virtual void Write_Frame(const std::string &frame_dir) const = 0;
	virtual void Save_Frame(const std::string &frame_dir) const = 0;
	virtual void Load_Frame(const std::string &frame_dir) = 0;

	virtual void Advance(const double dt) = 0;
};
