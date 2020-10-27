#include "GlApp.h"

#include <iostream>

#include <cstdlib>

namespace PhysX {

GlApp *GlApp::_this = nullptr;

GlApp::GlApp(const int width, const int height, const std::string &title) :
	_width(width),
	_height(height),
	_title(title)
{
	_this = this;
	// Initialize GLFW.
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, 4); // 4x MSAA
	// Create window.
	_window = glfwCreateWindow(_width, _height, _title.c_str(), NULL, NULL);
	if (!_window) {
		std::cerr << "Error: [GLApp] Failed to create GLFW window." << std::endl;
		std::exit(-1);
	}
	glfwMakeContextCurrent(_window);
	// Initialize GLAD.
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		std::cerr << "Error: [GLApp] Failed to initialize GLAD." << std::endl;
		std::exit(-1);
	}
}

void GlApp::run()
{
	initialize();
	while (!glfwWindowShouldClose(_window)) {
		processInput();

		update();
		clearBuffers();
		drawRenderItems();

		// Swap buffers and poll IO events.
		glfwSwapBuffers(_window);
		glfwPollEvents();
	}
}

void GlApp::initialize()
{
	setCallbacks();
	initPrograms();
	buildRenderItems();
}

void GlApp::setCallbacks() const
{
	glfwSetFramebufferSizeCallback(_window, onResize);
	glfwSetKeyCallback(_window, onKeyInput);
}

void GlApp::initPrograms()
{
	_programs["identity"] = std::make_unique<GlProgram>(_identityVsCode, _identityFsCode);
	_programs["flat"] = std::make_unique<GlProgram>(_flatVsCode, _flatFsCode);
}

void GlApp::clearBuffers() const
{
	glClearColor(_bgColor[0], _bgColor[1], _bgColor[2], 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

void GlApp::drawRenderItems() const
{
	for (auto &layer : _ritemLayers) {
		for (auto ritem : layer) {
			if (ritem->isVisible()) ritem->draw();
		}
	}
}

void GlApp::buildRenderItems()
{
	_ritems.push_back(std::make_unique<GlRenderTest>(_programs["identity"].get()));
	_ritemLayers[size_t(RenderLayer::Opaque)].push_back(_ritems.back().get());
}

void GlApp::processInput()
{
}

void GlApp::update()
{
	updateMsaaState();
	updateWireframeState();
}

void GlApp::onResize(GLFWwindow *window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void GlApp::onKeyInput(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS) {
		switch (key) {
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(window, true);
			break;
		case GLFW_KEY_SPACE:
			glfwSetWindowSize(window, _this->_width, _this->_height);
			break;
		case GLFW_KEY_F2:
			_this->_enableMsaa = !_this->_enableMsaa;
			break;
		case GLFW_KEY_F3:
			_this->_enableWireframe = !_this->_enableWireframe;
			break;
		}
	}
}

}
