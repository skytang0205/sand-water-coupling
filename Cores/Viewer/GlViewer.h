#pragma once
#pragma warning (disable : 4996)

#include "Graphics/GlApp.h"

#include <yaml-cpp/yaml.h>

#include <fmt/core.h>

#include <iostream>

#include <cstdio>

namespace PhysX {

class GlViewer : public GlApp
{
protected:

	const std::string &_outputDir;
	YAML::Node _root;

	uint _endFrame;
	uint _frameRate;
	uint _currentFrame = 0;

public:

	GlViewer(const std::string &outputDir, const uint frameRate);

	GlViewer() = delete;
	GlViewer(const GlViewer &rhs) = delete;
	GlViewer &operator=(const GlViewer &rhs) = delete;
	virtual ~GlViewer() = default;

protected:

	virtual void updateText() override;

	template <typename Type>
	Type loadItem(const std::string &name)
	{
		if (_root[name]) return _root[name].as<Type>();
		else {
			std::cerr << fmt::format("Error: [GlViewer] cannot find {} in description.yaml.", name) << std::endl;
			std::exit(-1);
		}
	}
};

}
