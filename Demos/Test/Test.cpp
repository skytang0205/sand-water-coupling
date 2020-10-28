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

inline void testEigen()
{
	using namespace PhysX;
	using namespace std::chrono;
	using namespace std::chrono_literals;
	Vector3f a = Vector3f::Zero();
	Vector3f b = Vector3f::Ones();
	Vector4f c = Vector4f::Zero();
	Vector4f d = Vector4f::Ones();
	const int T = 1000000000;
	auto t1 = steady_clock::now();
	for (int i = 0; i < T; i++) {
		a += (b * T).normalized();
	}
	auto t2 = steady_clock::now();
	for (int i = 0; i < T; i++) {
		c += (d * i).normalized();
	}
	auto t3 = steady_clock::now();
	std::cout << duration<double>(t2 - t1).count() << std::endl << duration<double>(t3 - t2).count() << std::endl;
}

int main(int argc, char *argv[])
{
	testEigen();
	return 0;
}
