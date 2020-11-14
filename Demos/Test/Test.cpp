#include "Physics/GridBasedFluid.h"
#include "Physics/Simulator.h"
#include "Utilities/ArgsParser.h"
#include "Utilities/Types.h"

#include <fmt/core.h>

#include <iostream>
#include <chrono>

inline void testArgsParser(int argc, char *argv[])
{
	PhysX::ArgsParser argsParser;
	argsParser.addArgument<std::string>("output", 'o', "the output directory", "output");
	argsParser.addArgument<int>("first", 'f', "the first frame", -1);
	argsParser.addArgument<int>("last", 'l', "the last frame", 200);
	argsParser.addArgument<double>("rate", 'r', "the frame rate (frames per second)", 50.0);
	argsParser.addArgument<int>("step", 's', "the step rate (steps per frame)", 1);
	argsParser.addArgument<std::string>("config", 'c', "the configuration file", "spring_mass_system.yaml");
	argsParser.parse(argc, argv);
	auto output = std::any_cast<std::string>(argsParser.getValueByName("output"));
	auto first = std::any_cast<int>(argsParser.getValueByName("first"));
	auto last = std::any_cast<int>(argsParser.getValueByName("last"));
	auto rate = std::any_cast<double>(argsParser.getValueByName("rate"));
	auto step = std::any_cast<int>(argsParser.getValueByName("step"));
	auto config = std::any_cast<std::string>(argsParser.getValueByName("config"));

	std::cout << fmt::format("{}\n{}\n{}\n{}\n{}\n{}\n", output, first, last, rate, step, config);

}

int main(int argc, char *argv[])
{
	using namespace PhysX;
	StaggeredGrid<2> grid(0.1, Vector2i(100, 100));
	auto fluid = std::make_unique<GridBasedFluid<2>>(grid);
	auto simulator = std::make_unique<Simulator>("output", 0, 50, 50, 1, fluid.get());
	simulator->Simulate();
	return 0;
}
