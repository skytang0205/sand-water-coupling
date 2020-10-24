#include "GlApp.h"

#include <memory>

int main()
{
	auto glApp = std::make_unique<PhysX::GlApp>("PhysX Viewer");
	glApp->run();
	return 0;
}
