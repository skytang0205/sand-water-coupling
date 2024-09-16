#include "Simulator.h"

#include "Utilities/Yaml.h"

#include <fmt/core.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace PhysX {

void Simulator::Simulate()
{
	using namespace std::chrono;
	const auto initialTime = steady_clock::now();

	if (_beginFrame == 0) {
		{ // Initialize simulation.
			std::cout << "[Initialize] ";
			_simulation->initialize();
			const auto endTime = steady_clock::now();
			std::cout << fmt::format(
				"    ... {:>7.2f}s used",
				duration<double>(endTime - initialTime).count())
				<< std::endl;
		}

		createOutputDirectory();
		writeAndSaveToFrameDirectory(0, true);

		{ // Write description.
			YAML::Node root;
			root["dimension"] = _simulation->dimension();
			_simulation->writeDescription(root);
			std::ofstream fout(_outputDir + "/description.yaml");
			fout << root;
		}

		_beginFrame = 1;
	}
	else {
		loadFromFrameDirectory(_beginFrame - 1);
	}

	// Initialize timing.
	const auto beginTime = steady_clock::now();
	std::cout << fmt::format("#  Time: {:>9.2f}s\n", duration<double>(beginTime - initialTime).count()) << std::endl;
	auto lastTime = beginTime;
	// Run the simulator.
	for (uint frame = _beginFrame; frame < _endFrame; frame++) {
		std::cout << fmt::format("** Simulate Frame {}...", frame) << std::endl;
		// Simulate.
		_simulation->setTime(real(frame - 1) / _frameRate);
		advanceTimeBySteps(real(1) / _frameRate);
		// Write and save files for frame.
		writeAndSaveToFrameDirectory(frame);
		// Output timing.
		const auto currentTime = steady_clock::now();
		std::cout << fmt::format(
			"#  Time: {:>9.2f}s / {:>9.2f}s\n"
			"   Prediction:        {:>9.2f}s\n",
			duration<double>(currentTime - lastTime).count(),
			duration<double>(currentTime - initialTime).count(),
			duration<double>(currentTime - beginTime).count() / (frame - _beginFrame + 1.0) * (_endFrame - _beginFrame)
				+ duration<double>(beginTime - initialTime).count()
			) << std::endl;
		lastTime = currentTime;
	}
}

void Simulator::createOutputDirectory() const
{
	std::cout << "** Create output directory..." << std::endl;
	std::filesystem::create_directories(_outputDir);
}

void Simulator::writeAndSaveToFrameDirectory(const uint frame, const bool staticDraw) const
{
	// Create, write and save to the frame directory.
	std::cout << fmt::format("** Write output files for frame {}...", frame) << std::endl;
	const std::string frameDir = fmt::format("{}/{}", _outputDir, frame);
	std::filesystem::create_directory(frameDir);
	_simulation->writeFrame(frameDir, staticDraw);
	_simulation->saveFrame(frameDir);
	// Write the last frame.
	std::ofstream fout(_outputDir + "/end_frame.txt");
	fout << frame + 1 << std::endl;
}

void Simulator::loadFromFrameDirectory(const uint frame)
{
	std::cout << fmt::format("** Load output files for frame {}...", frame) << std::endl;
	const std::string frameDir = fmt::format("{}/{}", _outputDir, frame);
	if (!std::filesystem::exists(frameDir)) {
		std::cerr << fmt::format("Error: [Simulator] No output directory for frame {}.", frame) << std::endl;
		std::exit(1);
	}
	_simulation->loadFrame(frameDir);
}

void Simulator::advanceTimeBySteps(const real targetTime)
{
	using namespace std::chrono;
	real time = 0.0;
	bool done = (time >= targetTime);
	while (!done) {
		const auto beginTime = steady_clock::now();

		real dt = _simulation->getTimeStep(_frameRate, _stepRate);
		if (time + dt >= targetTime) {
			dt = targetTime - time;
			done = true;
		}
		else if (time + 2 * dt >= targetTime) {
			dt = (targetTime - time) / 2;
		}

		// std::cout << fmt::format("[{:>6.0f} SPF] ", targetTime / dt);
		_simulation->advance(dt);
		time += dt;
		const auto endTime = steady_clock::now();
		// std::cout << fmt::format(
		// 	"    ... {:>7.2f}s used",
		// 	duration<double>(endTime - beginTime).count())
		// 	<< std::endl;
	}
}

}
