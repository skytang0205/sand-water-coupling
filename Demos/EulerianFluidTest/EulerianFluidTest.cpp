#include "EulerianFluidBuilder.h"

#include "Physics/Simulator.h"

int main(int argc, char *argv[])
{
	using namespace PhysX;
	auto fluid = EulerianFluidBuilder::build<2>(200);
	auto simulator = std::make_unique<Simulator>("output", 0, 100, 50, 1, fluid.get());
	simulator->Simulate();
	return 0;
}
