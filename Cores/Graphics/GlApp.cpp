#include "GlApp.h"

#include <fmt/core.h>

#include <iostream>

#include <cstdlib>

namespace PhysX {

GlApp *GlApp::_this = nullptr;

GlApp::GlApp(const int width, const int height, const std::string &title) :
	_savedWidth(width),
	_savedHeight(height),
	_savedTitle(title)
{
	_this = this;
	// Initialize GLFW.
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, 4); // 4x MSAA
	// Create window.
	_window = glfwCreateWindow(_savedWidth, _savedHeight, _savedTitle.c_str(), NULL, NULL);
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
	resize(_savedWidth, _savedHeight);
	initGlStates();
	setCallbacks();
	initPrograms();
	buildRenderItems();
}

void GlApp::resize(const int width, const int height)
{
	glViewport(0, 0, width, height);
	_orbitCamera.setAspectRatio(float(width) / height);
}

void GlApp::initGlStates() const
{
	glEnable(GL_FRAMEBUFFER_SRGB); // gamma correction
}

void GlApp::setCallbacks() const
{
	glfwSetFramebufferSizeCallback(_window, framebufferSizeCallback);
	glfwSetKeyCallback(_window, keyCallback);
	glfwSetMouseButtonCallback(_window, mouseButtonCallback);
	glfwSetCursorPosCallback(_window, cursorPosCallback);
	glfwSetScrollCallback(_window, scrollCallback);
}

void GlApp::initPrograms()
{
	_programs["identity"] = std::make_unique<GlProgram>(_kIdentityVsCode, _kIdentityFsCode);
	_programs["flat"] = std::make_unique<GlProgram>(_kFlatVsCode, _kFlatFsCode);
}

void GlApp::buildRenderItems()
{
	_ritems.push_back(std::make_unique<GlRenderTest>(_programs["flat"].get()));
	_ritemLayers[size_t(RenderLayer::Opaque)].push_back(_ritems.back().get());
}

void GlApp::processInput()
{
}

void GlApp::update()
{
	updateMsaaState();
	updateWireframeState();
	_orbitCamera.update();
	updateUniforms();
	updateFrameRate();
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

void GlApp::updateFrameRate()
{
	static uint ticks = 0;
	static auto lastTime = glfwGetTime();

	ticks++;
	auto currentTime = glfwGetTime();
	auto elapsed = currentTime - lastTime;
	if (elapsed > 1) {
		std::cout << fmt::format("FPS = {:.0f}", ticks / elapsed) << std::endl;
		lastTime = currentTime;
		ticks = 0;
	}
}

void GlApp::updateUniforms()
{
	_programs["flat"]->setUniform("uWorld", _orbitCamera.projView());
}

void GlApp::framebufferSizeCallback(GLFWwindow *window, int width, int height)
{
	_this->resize(width, height);
}

void GlApp::keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS) {
		switch (key) {
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(window, true);
			break;
		case GLFW_KEY_SPACE:
			glfwSetWindowSize(window, _this->_savedWidth, _this->_savedHeight);
			_this->_orbitCamera.reset();
			break;
		case GLFW_KEY_F2:
			_this->_enableMsaa = (_this->_enableMsaa ^ 1) | 2;
			break;
		case GLFW_KEY_F3:
			_this->_enableWireframe = (_this->_enableWireframe ^ 1) | 2;
			break;
		}
	}
}

void GlApp::mouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
{
	if (action == GLFW_PRESS) {
		glfwGetCursorPos(window, &_this->_lastMousePos.x(), &_this->_lastMousePos.y());
	}
}

void GlApp::cursorPosCallback(GLFWwindow *window, double xpos, double ypos)
{
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
		_this->_orbitCamera.rotate(float(xpos - _this->_lastMousePos.x()), float(ypos - _this->_lastMousePos.y()));
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
		_this->_orbitCamera.translate(float(xpos - _this->_lastMousePos.x()), float(ypos - _this->_lastMousePos.y()));
	}
	_this->_lastMousePos << xpos, ypos;
}

void GlApp::scrollCallback(GLFWwindow *window, double xoffset, double yoffset)
{
	_this->_orbitCamera.scale(float(yoffset));
}

}
