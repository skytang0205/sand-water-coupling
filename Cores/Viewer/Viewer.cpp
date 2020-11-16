#include "GlViewer.h"

#include "Utilities/ArgsParser.h"

inline std::unique_ptr<PhysX::ArgsParser> BuildArgsParser()
{
	auto parser = std::make_unique<PhysX::ArgsParser>();
	parser->addArgument<std::string>("output", 'o', "the output directory", "output");
	parser->addArgument<int>("rate", 'r', "the frame rate (frames per second)", 50);
	return parser;
}

int main(int argc, char *argv[])
{
	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

	const auto output = std::any_cast<std::string>(parser->getValueByName("output"));
	const auto rate = std::any_cast<int>(parser->getValueByName("rate"));

	auto glApp = std::make_unique<PhysX::GlViewer>(output, rate);
	glApp->run();
	return 0;
}
