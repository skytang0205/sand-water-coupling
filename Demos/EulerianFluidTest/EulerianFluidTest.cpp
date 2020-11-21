#include "EulerianFluidBuilder.h"

#include "Physics/Simulator.h"
#include "Utilities/ArgsParser.h"

using namespace PhysX;

inline std::unique_ptr<ArgsParser> BuildArgsParser()
{
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("output", 'o', "the output directory", "output");
	parser->addArgument<int>("test", 't', "the test case index", 0);
	parser->addArgument<uint>("begin", 'b', "the begin frame (including)", 0);
	parser->addArgument<uint>("end", 'e', "the end frame (excluding)", 200);
	parser->addArgument<uint>("rate", 'r', "the frame rate (frames per second)", 50);
	parser->addArgument<real>("cfl", 'c', "the CFL number", real(1));
	parser->addArgument<int>("scale", 's', "the scale of grid", 200);
	return parser;
}

int main(int argc, char *argv[])
{
	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

	const auto output = std::any_cast<std::string>(parser->getValueByName("output"));
	const auto test = std::any_cast<int>(parser->getValueByName("test"));
	const auto begin = std::any_cast<uint>(parser->getValueByName("begin"));
	const auto end = std::any_cast<uint>(parser->getValueByName("end"));
	const auto rate = std::any_cast<uint>(parser->getValueByName("rate"));
	const auto cfl = std::any_cast<real>(parser->getValueByName("cfl"));
	const auto scale = std::any_cast<int>(parser->getValueByName("scale"));

	auto fluid = PhysX::EulerianFluidBuilder::build<2>(scale);
	auto simulator = std::make_unique<PhysX::Simulator>(output, begin, end, rate, cfl, fluid.get());
	simulator->Simulate();

	return 0;
}
