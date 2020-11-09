#pragma once

#include "Utilities/Types.h"

#include <fstream>
#include <string>

namespace PhysX {

class Simulation
{
protected:

	real _time = 0;

public:

	Simulation() = default;
	Simulation(const Simulation &rhs) = delete;
	Simulation &operator=(const Simulation &rhs) = delete;
	virtual ~Simulation() = default;

	void setTime(const real time) { _time = time; }
	virtual real getTimeStep(const uint frameRate, const real stepRate) const = 0;

	virtual void writeDescription(std::ofstream &output) const = 0;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const = 0;
	virtual void saveFrame(const std::string &frameDir) const = 0;
	virtual void loadFrame(const std::string &frameDir) = 0;

	virtual void advance(const real dt) = 0;
};

}
