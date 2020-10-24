#pragma once

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <string>

namespace PhysX {

class GlApp
{
protected:

	GLFWwindow *_window;

public:

	GlApp(const std::string &caption);

	GlApp() = delete;
	GlApp(const GlApp &rhs) = delete;
	GlApp &operator=(const GlApp &rhs) = delete;
	virtual ~GlApp() { glfwTerminate(); }

	void Run();

	void processInput();

private:

	static void onResize(GLFWwindow *window, int width, int height);
};

}
