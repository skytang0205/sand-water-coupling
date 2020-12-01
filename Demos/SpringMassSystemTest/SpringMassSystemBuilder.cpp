#include "SpringMassSystemBuilder.h"

namespace PhysX {

void SpringMassSystemBuilder::reportError(const std::string &msg)
{
	std::cerr << fmt::format("Error: [EulerianFluidBuilder] encountered {}.\n{}", msg) << std::endl;
	std::exit(-1);
}

}
