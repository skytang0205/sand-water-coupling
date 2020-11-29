#include "ParticleInCellLiquidBuilder.h"

#include "Physics/Simulator.h"
#include "Utilities/ArgsParser.h"

using namespace PhysX;

inline std::unique_ptr<ArgsParser> BuildArgsParser()
{
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("output", 'o', "the output directory", "output");
	parser->addArgument<int>("dim", 'd', "the dimension of the simulation", 2);
	parser->addArgument<int>("test", 't', "the test case index", 0);
	parser->addArgument<uint>("begin", 'b', "the begin frame (including)", 0);
	parser->addArgument<uint>("end", 'e', "the end frame (excluding)", 200);
	parser->addArgument<uint>("rate", 'r', "the frame rate (frames per second)", 50);
	parser->addArgument<real>("cfl", 'c', "the CFL number", real(1));
	parser->addArgument<int>("scale", 's', "the scale of grid", -1);
	parser->addArgument<int>("markers", 'm', "the number of markers per sub-cell", 2);
	parser->addArgument<real>("alpha", 'a', "the proportion of PIC algorithm versus FLIP", 0);
	return parser;
}

int main(int argc, char *argv[])
{
#ifdef _OPENMP
	omp_set_num_threads(std::max(omp_get_num_procs() / 3, 1));
#endif

	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

	const auto output = std::any_cast<std::string>(parser->getValueByName("output"));
	const auto dim = std::any_cast<int>(parser->getValueByName("dim"));
	const auto test = std::any_cast<int>(parser->getValueByName("test"));
	const auto begin = std::any_cast<uint>(parser->getValueByName("begin"));
	const auto end = std::any_cast<uint>(parser->getValueByName("end"));
	const auto rate = std::any_cast<uint>(parser->getValueByName("rate"));
	const auto cfl = std::any_cast<real>(parser->getValueByName("cfl"));
	const auto scale = std::any_cast<int>(parser->getValueByName("scale"));
	const auto markers = std::any_cast<int>(parser->getValueByName("markers"));
	const auto alpha = std::any_cast<real>(parser->getValueByName("alpha"));

	std::unique_ptr<Simulation> simulation;
	if (dim == 2)
		simulation = ParticleInCellLiquidBuilder::build<2>(scale, test, markers, alpha);
	else if (dim == 3)
		simulation = ParticleInCellLiquidBuilder::build<3>(scale, test, markers, alpha);
	else {
		std::cerr << "Error: [main] encountered invalid dimension." << std::endl;
		std::exit(-1);
	}
	auto simulator = std::make_unique<Simulator>(output, begin, end, rate, cfl, simulation.get());
	simulator->Simulate();

	return 0;
}
