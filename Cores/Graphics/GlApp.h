#pragma once

#include "GlRenderItem.h"
#include "Types.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace PhysX {

class GlApp
{
public:

	enum class RenderLayer : uint { Opaque, Text, Count };

protected:

	GLFWwindow *_window;

	const std::string &_title;
	Vector3f _bgColor = Vector3f(176, 196, 222) / 255;

	std::unordered_map<std::string, std::unique_ptr<GlProgram>> _programs;

	std::vector<std::unique_ptr<GlRenderItem>> _ritems;
	std::vector<GlRenderItem *> _ritemLayer[size_t(RenderLayer::Count)];

public:

	GlApp(const int width, const int height, const std::string &caption);

	GlApp() = delete;
	GlApp(const GlApp &rhs) = delete;
	GlApp &operator=(const GlApp &rhs) = delete;
	virtual ~GlApp() { glfwTerminate(); }

	void run();

	virtual void buildRenderItems();
	virtual void processInput();
	virtual void update();
	virtual void clearBuffers() const;
	virtual void drawRenderItems() const; // Render two triangles (only for test).

private:

	static void onResize(GLFWwindow *window, int width, int height);
};

}
