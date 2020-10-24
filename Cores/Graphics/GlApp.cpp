#include "GLApp.h"

#include <iostream>

#include <cstdlib>

namespace PhysX {

GlApp::GlApp(const std::string &caption)
{
	// Initialize GLFW.
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	// Initialize window.
	_window = glfwCreateWindow(1024, 768, caption.c_str(), NULL, NULL);
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
}

void GlApp::run()
{
	while (!glfwWindowShouldClose(_window)) {
		processInput();

		// Render.
		glClearColor(176 / 255.0, 196 / 255.0, 222 / 255.0, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		// Swap buffers and poll IO events.
		glfwSwapBuffers(_window);
		glfwPollEvents();
	}
}

void GlApp::processInput()
{
	if (glfwGetKey(_window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(_window, true);
}

void GlApp::onResize(GLFWwindow *window, int width, int height)
{
	glViewport(0, 0, width, height);
}

}
