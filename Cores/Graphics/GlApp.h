#pragma once

#include "Graphics/GlRenderItem.h"
#include "Utilities/Types.h"

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

private:

	static GlApp *_this;

	const GLchar *_identityVsCode =
#include "GlIdentityShader.vert"
		;
	const GLchar *_identityFsCode =
#include "GlIdentityShader.frag"
		;
	const GLchar *_flatVsCode =
#include "GlFlatShader.vert"
		;
	const GLchar *_flatFsCode =
#include "GlFlatShader.frag"
		;

	bool _enableMsaa = false;
	bool _enableWireframe = false;

protected:

	GLFWwindow *_window;

	const std::string &_title;
	Vector3f _bgColor = Vector3f(176, 196, 222) / 255;

	std::unordered_map<std::string, std::unique_ptr<GlProgram>> _programs;

	std::vector<std::unique_ptr<GlRenderItem>> _ritems;
	std::vector<GlRenderItem *> _ritemLayers[size_t(RenderLayer::Count)];

public:

	GlApp(const int width, const int height, const std::string &caption);

	GlApp() = delete;
	GlApp(const GlApp &rhs) = delete;
	GlApp &operator=(const GlApp &rhs) = delete;
	virtual ~GlApp() { glfwTerminate(); }

	void run();

protected:

	virtual void initPrograms();
	virtual void buildRenderItems();
	virtual void processInput();
	virtual void update();
	virtual void clearBuffers() const;
	virtual void drawRenderItems() const;

private:

	void updateMsaaState() { _enableMsaa ? glEnable(GL_MULTISAMPLE) : glDisable(GL_MULTISAMPLE); }
	void updateWireframeState() { _enableWireframe ? glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) : glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); }

	static void onResize(GLFWwindow *window, int width, int height);
	static void onKeyInput(GLFWwindow *window, int key, int scancode, int action, int mods);
};

}
