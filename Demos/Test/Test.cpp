#include "Utilities/ArgsParser.h"
#include "Utilities/Types.h"

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

struct PassConstants
{
	PhysX::Matrix4f projView;			// 0
	PhysX::Vector3f viewPos;			// 64
	float pad0;
	PhysX::Vector3f ambientStrength;	// 80
	float pad1;
	PhysX::Vector3f lightStrength;		// 96
	float pad2;
	PhysX::Vector3f lightDir;			// 112
	float pad3;
	float totalTime;			// 128
	float deltaTime;			// 132
} test;

inline void testEigen()
{
	std::cout << (char *)(&test.viewPos) - (char *)(&test) << std::endl;
	std::cout << (char *)(&test.ambientStrength) - (char *)(&test) << std::endl;
	std::cout << (char *)(&test.lightStrength) - (char *)(&test) << std::endl;
	std::cout << (char *)(&test.lightDir) - (char *)(&test) << std::endl;
	std::cout << (char *)(&test.totalTime) - (char *)(&test) << std::endl;
	std::cout << (char *)(&test.deltaTime) - (char *)(&test) << std::endl;
}

int main(int argc, char *argv[])
{
	testArgsParser(argc, argv);
	testEigen();
	return 0;
}
