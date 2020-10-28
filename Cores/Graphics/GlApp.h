#pragma once

#include "Graphics/GlOrbitCamera.h"
#include "Graphics/GlRenderItem.h"
#include "Utilities/Types.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <memory>
#include <numbers>
#include <string>
#include <unordered_map>
#include <vector>

namespace PhysX {

class GlApp
{
protected:

	enum class RenderLayer : uint { Opaque, Text, Count };

	static GlApp *_this;

	static constexpr GLchar _kIdentityVsCode[] =
#include "GlIdentityShader.vert"
		;
	static constexpr GLchar _kIdentityFsCode[] =
#include "GlIdentityShader.frag"
		;
	static constexpr GLchar _kFlatVsCode[] =
#include "GlFlatShader.vert"
		;
	static constexpr GLchar _kFlatFsCode[] =
#include "GlFlatShader.frag"
		;

	uchar _enableMsaa = 2;
	uchar _enableWireframe = 2;

	GLFWwindow *_window;

	const int _savedWidth;
	const int _savedHeight;
	const std::string &_savedTitle;

	Vector3f _bgColor = Vector3f(176, 196, 222) / 255;
	Vector2d _lastMousePos = Vector2d::Zero();

	GlOrbitCamera _orbitCamera = GlOrbitCamera(
		0.25f * float(std::numbers::pi), // fovy
		1.0f, // aspect
		1.0f, // zNear
		1000.0f, // zFar
		10.0f, // radius
		0.0, // phi
		0.5f * float(std::numbers::pi), // theta
		Vector3f::Zero() // target
		);

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

	virtual void resize(const int width, const int height);
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

	void updateUniforms();

	void updateFrameRate();

	static void framebufferSizeCallback(GLFWwindow *window, int width, int height);
	static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
	static void mouseButtonCallback(GLFWwindow *window, int button, int action, int mods);
	static void cursorPosCallback(GLFWwindow *window, double xpos, double ypos);
	static void scrollCallback(GLFWwindow *window, double xoffset, double yoffset);
};

}
