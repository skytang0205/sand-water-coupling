#include "GlViewer.h"

#include <algorithm>
#include <fstream>

namespace PhysX {

GlViewer::GlViewer(const std::string &outputDir, const uint frameRate) :
	GlApp(1024, 768, "PhysX Viewer - " + outputDir),
	_outputDir(outputDir),
	_frameRate(frameRate)
{
	_this = this;
	{ // Load desciption.yaml.
		std::ifstream fin(_outputDir + "/description.yaml");
		if (!fin) {
			std::cerr << fmt::format("Error: [GlViewer] failed to load {}/description.yaml.", _outputDir) << std::endl;
			std::exit(-1);
		}
		_root = YAML::Load(fin);
		if (_root["dimension"]) _dim = _root["dimension"].as<int>();
		else {
			std::cerr << "Error: [GlViewer] cannot find dimension in description.yaml." << std::endl;
			std::exit(-1);
		}
	}
	{ // load end_frame.txt
		std::ifstream fin(_outputDir + "/end_frame.txt");
		if (!fin) {
			std::cerr << fmt::format("Error: [GlViewer] failed to load {}/end_frame.txt.", _outputDir) << std::endl;
			std::exit(-1);
		}
		fin >> _endFrame;
	}
}

void GlViewer::setCallbacks() const
{
	GlApp::setCallbacks();
	glfwSetKeyCallback(_window, keyCallback);
}

void GlViewer::buildRenderItems()
{
	GlApp::buildRenderItems();

	if (!_root["objects"] || !_root["objects"].IsSequence()) {
		std::cerr << fmt::format("Error: [GlViewer] cannot find objects sequence in description.yaml.") << std::endl;
		std::exit(-1);
	}

	for (const auto &node : _root["objects"]) {
		auto _simulated = std::make_unique<GlSimulated>(_programs["default"].get(), _outputDir, _endFrame, _dim, node);

		// Push the item into vectors.
		_ritemLayers[uint(_simulated->isTransparent() ? RenderLayer::Transparency : RenderLayer::Opaque)].push_back(_simulated.get());
		_simulatedObjects.push_back(_simulated.get());
		_ritems.push_back(std::move(_simulated));
	}
}

void GlViewer::update(const double dt)
{
	if (_playing) {
		_currentFrame += 1.0 / (dt * _frameRate);
		if (_currentFrame >= _endFrame - 1) {
			_currentFrame = _endFrame - 1;
			_playing = false;
		}
	}
	for (auto simulated : _simulatedObjects)
		simulated->setCurrentFrame(uint(_currentFrame));
	GlApp::update(dt);
}

void GlViewer::updateText()
{
	GlApp::updateText();
	_text->set(
		fmt::format("Frame: {:>4}", uint(_currentFrame)),
		Vector2f(1.0f - 10.24f / _width, 1.0f - (24.0f + 7.68f) / _height),
		Vector2f(0.75f, 0.75f),
		Vector4f(0, 0, 0, 1),
		GlText::Alignment::Right
	);
}

void GlViewer::keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS) {
		switch (key) {
		case GLFW_KEY_P:
			_this->_playing = !_this->_playing;
			break;
		case GLFW_KEY_R:
			_this->_currentFrame = 0;
			break;
		case GLFW_KEY_LEFT_BRACKET:
			if (!_this->_playing) {
				_this->_currentFrame = std::max(_this->_currentFrame - 1.0, 0.0);
			}
			break;
		case GLFW_KEY_RIGHT_BRACKET:
			if (!_this->_playing) {
				_this->_currentFrame = std::min(_this->_currentFrame + 1.0, _this->_endFrame - 1.0);
			}
			break;
		case GLFW_KEY_1:
		case GLFW_KEY_2:
		case GLFW_KEY_3:
		case GLFW_KEY_4:
		case GLFW_KEY_5:
		case GLFW_KEY_6:
		case GLFW_KEY_7:
		case GLFW_KEY_8:
		case GLFW_KEY_9:
			if (key - size_t('1') < _this->_simulatedObjects.size())
				_this->_simulatedObjects[key - size_t('1')]->flipVisibility();
			break;
		default:
			GlApp::keyCallback(window, key, scancode, action, mods);
			break;
		}
	}
	else if (action == GLFW_REPEAT) {
		switch (key) {
		case GLFW_KEY_LEFT_BRACKET:
			if (!_this->_playing) {
				_this->_currentFrame = std::max(_this->_currentFrame - 1.0, 0.0);
			}
			break;
		case GLFW_KEY_RIGHT_BRACKET:
			if (!_this->_playing) {
				_this->_currentFrame = std::min(_this->_currentFrame + 1.0, _this->_endFrame - 1.0);
			}
			break;
		default:
			GlApp::keyCallback(window, key, scancode, action, mods);
			break;
		}
	}
}

}

