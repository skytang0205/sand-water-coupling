#include "LevelSetLiquidBuilder.h"

namespace PhysX {

void LevelSetLiquidBuilder::reportError(const std::string &msg)
{
	std::cerr << fmt::format("Error: [LevelSetLiquidBuilder] encountered {}.\n{}", msg) << std::endl;
	std::exit(-1);
}

}
