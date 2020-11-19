#include "LevelSetLiquidBuilder.h"

#include "Physics/Simulator.h"
#include "Utilities/ArgsParser.h"

inline std::unique_ptr<PhysX::ArgsParser> BuildArgsParser()
{
	auto parser = std::make_unique<PhysX::ArgsParser>();
	parser->addArgument<std::string>("output", 'o', "the output directory", "output");
	parser->addArgument<int>("begin", 'b', "the begin frame (including)", 0);
	parser->addArgument<int>("end", 'e', "the end frame (excluding)", 200);
	parser->addArgument<int>("rate", 'r', "the frame rate (frames per second)", 50);
	parser->addArgument<PhysX::real>("cfl", 'c', "the CFL number", PhysX::real(1));
	parser->addArgument<int>("scale", 's', "the scale of grid", 200);
	return parser;
}

int main(int argc, char *argv[])
{
	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

	const auto output = std::any_cast<std::string>(parser->getValueByName("output"));
	const auto begin = std::any_cast<int>(parser->getValueByName("begin"));
	const auto end = std::any_cast<int>(parser->getValueByName("end"));
	const auto rate = std::any_cast<int>(parser->getValueByName("rate"));
	const auto cfl = std::any_cast<PhysX::real>(parser->getValueByName("cfl"));
	const auto scale = std::any_cast<int>(parser->getValueByName("scale"));

	auto liquid = PhysX::LevelSetLiquidBuilder::build<2>(scale);
	auto simulator = std::make_unique<PhysX::Simulator>(output, begin, end, rate, cfl, liquid.get());
	simulator->Simulate();

	auto box = PhysX::ImplicitBox<2>(PhysX::Vector2r(1, 2), PhysX::Vector2r(4, 3));
	std::cout << box.closestPosition(PhysX::Vector2r(0, 0)).transpose() << std::endl;
	std::cout << box.closestPosition(PhysX::Vector2r(5, 5)).transpose() << std::endl;
	std::cout << box.closestPosition(PhysX::Vector2r(-2, 10)).transpose() << std::endl;
	std::cout << box.closestPosition(PhysX::Vector2r(-6, -6)).transpose() << std::endl;

	return 0;
}
