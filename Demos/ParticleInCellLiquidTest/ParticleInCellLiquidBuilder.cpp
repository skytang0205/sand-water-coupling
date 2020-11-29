#include "ParticleInCellLiquidBuilder.h"

namespace PhysX {

void ParticleInCellLiquidBuilder::reportError(const std::string &msg)
{
	std::cerr << fmt::format("Error: [ParticleInCellLiquidBuilder] encountered {}.\n{}", msg) << std::endl;
	std::exit(-1);
}

}
