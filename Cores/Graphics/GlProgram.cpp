#include "GlProgram.h"

#include <fmt/core.h>

#include <iostream>

#include <cstdlib>

namespace PhysX {

 GlProgram::GlProgram(const GLchar *vsCode, const GLchar *fsCode)
{
	// Compile VS.
	auto vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vsCode, nullptr);
	glCompileShader(vertexShader);
	checkCompileErrors(vertexShader, "vertex");
	// Compile FS.
	auto fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fsCode, nullptr);
	glCompileShader(fragmentShader);
	checkCompileErrors(fragmentShader, "fragment");
	// Link shaders to program.
	_program = glCreateProgram();
	glAttachShader(_program, vertexShader);
	glAttachShader(_program, fragmentShader);
	glLinkProgram(_program);
	checkCompileErrors(_program, "program");
	// Delete the shaders.
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
}

void GlProgram::checkCompileErrors(GLuint object, const std::string &type) const
{
	static const GLsizei kLogSize = 1024;
	GLchar *infoLog;

	GLint success;
	if (type != "program") {
		glGetShaderiv(object, GL_COMPILE_STATUS, &success);
		if (!success) {
			infoLog = new GLchar[kLogSize];
			glGetShaderInfoLog(object, kLogSize, nullptr, infoLog);
			std::cerr << fmt::format("Error: [GlShader] failed to compile {} shader.\n{}", type, infoLog) << std::endl;
		}
	}
	else {
		glGetProgramiv(object, GL_LINK_STATUS, &success);
		if (!success) {
			infoLog = new GLchar[kLogSize];
			glGetProgramInfoLog(object, kLogSize, nullptr, infoLog);
			std::cerr << fmt::format("Error: [GlShader] failed to link {}.\n{}", type, infoLog) << std::endl;
		}
	}
	std::exit(-1);
}

}

