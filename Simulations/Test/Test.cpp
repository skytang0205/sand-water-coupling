#include "ArgsParser.h"

#include <iostream>

int main()
{
	PhysX::ArgsParser argsParser;
	argsParser.addArgument<std::string>("output", 'o', "the output directory", std::string("output"));
	argsParser.addArgument<int>("first", 'f', "the first frame", -1);
	argsParser.addArgument<int>("last", 'l', "the last frame", 200);
	argsParser.addArgument<double>("rate", 'r', "the frame rate (frames per second)", 50.0);
	argsParser.addArgument<int>("step", 's', "the step rate (steps per frame)", 1);
	argsParser.addArgument<std::string>("config", 'c', "the configuration file", std::string("spring_mass_system.yaml"));
	argsParser.addArgument<int>("ddd");
	argsParser.parse("Test -j 3");
	return 0;
}
