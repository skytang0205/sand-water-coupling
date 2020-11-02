#include "Graphics/GeometryGenerator.h"
#include "Graphics/GlApp.h"

#include <memory>

int main()
{
	PhysX::GeometryGenerator::createBox(1.0, 2.0, 3.0);
	auto glApp = std::make_unique<PhysX::GlApp>(1024, 768, "PhysX Viewer");
	glApp->run();
	return 0;
}
