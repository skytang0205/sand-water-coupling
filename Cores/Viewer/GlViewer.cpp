#include "GlViewer.h"

#include <fstream>

namespace PhysX {

GlViewer::GlViewer(const std::string &outputDir, const uint frameRate) :
	GlApp(1024, 768, "PhysX Viewer - " + outputDir),
	_outputDir(outputDir),
	_frameRate(frameRate)
{
	{ // Load desciption.yaml.
		std::ifstream input(_outputDir + "/description.yaml");
		if (!input) {
			std::cerr << fmt::format("Error: [GlViewer] failed to load {}/description.yaml.", _outputDir) << std::endl;
			std::exit(-1);
		}
		_root = YAML::Load(input);
		_dim = loadItem<int>("dimension");
	}
	{ // load end_frame.txt
		std::ifstream input(_outputDir + "/end_frame.txt");
		if (!input) {
			std::cerr << fmt::format("Error: [GlViewer] failed to load {}/end_frame.txt.", _outputDir) << std::endl;
			std::exit(-1);
		}
		input >> _endFrame;
	}
}

void GlViewer::updateText()
{
	GlApp::updateText();
	_text->set(
		fmt::format("Frame: {:>4}", _currentFrame),
		Vector2f(1.0f - 10.24f / _width, 1.0f - (24.0f + 7.68f) / _height),
		Vector2f(0.75f, 0.75f),
		Vector4f(0, 0, 0, 1),
		GlText::Alignment::Right
	);
}

}

