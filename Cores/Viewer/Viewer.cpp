#include "GlViewer.h"

#include <memory>

int main()
{
	auto glApp = std::make_unique<PhysX::GlViewer>("output", 50);
	glApp->run();
	return 0;
}
