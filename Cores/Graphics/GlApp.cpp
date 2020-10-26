#include "GlApp.h"

#include <iostream>

#include <cstdlib>

namespace PhysX {

GlApp::GlApp(const int width, const int height, const std::string &title) : _title(title)
{
	// Initialize GLFW.
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	// Initialize window.
	_window = glfwCreateWindow(width, height, _title.c_str(), NULL, NULL);
	if (_window == nullptr) {
		std::cerr << "Error: [GLApp] Failed to create GLFW window." << std::endl;
		std::exit(-1);
	}
	glfwMakeContextCurrent(_window);
	glfwSetFramebufferSizeCallback(_window, onResize);
	// Initialize GLAD.
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		std::cerr << "Error: [GLApp] Failed to initialize GLAD." << std::endl;
		std::exit(-1);
	}

	_programs["identity"] = std::make_unique<GlProgram>(_identityVsCode, _identityFsCode);

	_ritems.push_back(std::make_unique<GlRenderTest>(_programs["identity"].get()));
	_ritemLayers[size_t(RenderLayer::Opaque)].push_back(_ritems.back().get());
}

void GlApp::run()
{
	buildRenderItems();
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

void GlApp::clearBuffers() const
{
	glClearColor(_bgColor[0], _bgColor[1], _bgColor[2], 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

void GlApp::drawRenderItems() const
{
	for (auto &layer : _ritemLayers) {
		for (auto ritem : layer) {
			ritem->draw();
		}
	}
}

void GlApp::buildRenderItems()
{ }

void GlApp::processInput()
{
	if (glfwGetKey(_window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(_window, true);
}

void GlApp::update()
{ }

void GlApp::onResize(GLFWwindow *window, int width, int height) 
{
	glViewport(0, 0, width, height);
}

}
