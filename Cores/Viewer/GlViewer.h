#pragma once
#pragma warning (disable : 4996)

#include "GlSimulated.h"

#include "Graphics/GlApp.h"
#include "Utilities/Yaml.h"

#include <fmt/core.h>

#include <iostream>

#include <cstdio>

namespace PhysX {

class GlViewer : public GlApp
{
protected:

	static inline GlViewer *_this = nullptr;

	const std::string &_outputDir;
	YAML::Node _root;

	uint _endFrame;
	uint _frameRate;

	bool _playing = false;
	double _currentFrame = 0;

	std::vector<GlSimulated *> _simulatedObjects;

public:

	GlViewer(const std::string &outputDir, const uint frameRate);

	GlViewer() = delete;
	GlViewer(const GlViewer &rhs) = delete;
	GlViewer &operator=(const GlViewer &rhs) = delete;
	virtual ~GlViewer() = default;

protected:

	virtual void setCallbacks() const override;
	virtual void buildRenderItems() override;

	virtual void update(const double dt) override;
	virtual void updateText() override;

	static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
};

}
