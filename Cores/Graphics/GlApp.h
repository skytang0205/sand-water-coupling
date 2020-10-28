#pragma once

#include "Graphics/GlCamera.h"
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
protected:

	enum class RenderLayer : uint { Opaque, Text, Count };

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


	uchar _enableMsaa = 2;
	uchar _enableWireframe = 2;

	GLFWwindow *_window;

	const int _width; // the default width
	const int _height; // the default height
	const std::string &_title;

	Vector3f _bgColor = Vector3f(176, 196, 222) / 255;
	Vector2d _lastMousePos = Vector2d::Zero();

	GlPolarCamera _polarCamera;

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

	void initialize();

	virtual void initGlStates() const;
	virtual void setCallbacks() const;
	virtual void initPrograms();
	virtual void buildRenderItems();

	virtual void processInput();
	virtual void update();
	virtual void clearBuffers() const;
	virtual void drawRenderItems() const;

	void updateMsaaState()
	{
		if (_enableMsaa & 2)
			_enableMsaa & 1 ? glEnable(GL_MULTISAMPLE) : glDisable(GL_MULTISAMPLE);
		_enableMsaa &= 1;
	}

	void updateWireframeState()
	{
		if (_enableWireframe & 2)
			_enableWireframe & 1 ? glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) : glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		_enableWireframe &= 1;
	}

	void updateFrameRate();

	static void framebufferSizeCallback(GLFWwindow *window, int width, int height);
	static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
	static void mouseButtonCallback(GLFWwindow *window, int button, int action, int mods);
	static void cursorPosCallback(GLFWwindow *window, double xpos, double ypos);
	static void scrollCallback(GLFWwindow *window, double xoffset, double yoffset);
};

}
