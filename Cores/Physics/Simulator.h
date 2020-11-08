#pragma once

#include "Physics/Simulation.h"
#include "Utilities/Types.h"

#include <string>

namespace PhysX {

class Simulator
{
protected:

	const std::string _outputDir;
	uint _beginFrame;
	uint _endFrame;
	const uint _frameRate; // number of frames per second
	const real _stepRate; // number of steps per frame

	Simulation *const _simulation;

public:

	Simulator(
		const std::string &outputDir,
		const uint beginFrame,
		const uint endFrame,
		const uint frameRate,
		const double stepRate,
		Simulation *const simulation
		) :
		_outputDir(outputDir),
		_beginFrame(beginFrame),
		_endFrame(endFrame),
		_frameRate(frameRate),
		_stepRate(stepRate),
		_simulation(simulation)
	{ }

	Simulator() = delete;
	Simulator(const Simulator &rhs) = delete;
	Simulator &operator=(const Simulator &rhs) = delete;
	virtual ~Simulator() = default;

	void Simulate();

protected:

	void createOutputDirectory() const;
	void writeAndSaveToFrameDirectory(const uint frame, const bool staticDraw) const;
	void loadFromFrameDirectory(const uint frame);
	void advanceTimeBySteps(const real targetTime);
};

}
