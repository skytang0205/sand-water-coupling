#include "GlApp.h"

#include <memory>

int main()
{
	auto glApp = std::make_unique<PhysX::GlApp>(1024, 768, "PhysX Viewer");
	glApp->run();
	return 0;
}
