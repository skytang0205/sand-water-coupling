#pragma once

#include "Types.h"

#include <glad/glad.h>

#include <string>

namespace PhysX {

class GlProgram
{
protected:

	GLuint _program;

public:

	GlProgram(const std::string &vsFileName, const std::string &fsFileName);

	GlProgram() = delete;
	GlProgram(const GlProgram &rhs) = delete;
	GlProgram &operator=(const GlProgram &rhs) = delete;
	virtual ~GlProgram() { }

	void use() const { glUseProgram(_program); }

	void setUniform(const std::string &name, const GLfloat value) const { glUniform1f(glGetUniformLocation(_program, name.c_str()), value); }
	void setUniform(const std::string &name, const Vector2f value) const { glUniform2f(glGetUniformLocation(_program, name.c_str()), value[0], value[1]); }
	void setUniform(const std::string &name, const Vector3f value) const { glUniform3f(glGetUniformLocation(_program, name.c_str()), value[0], value[1], value[2]); }
	void setUniform(const std::string &name, const Vector4f value) const { glUniform4f(glGetUniformLocation(_program, name.c_str()), value[0], value[1], value[2], value[3]); }

	void setUniform(const std::string &name, const GLint value) const { glUniform1i(glGetUniformLocation(_program, name.c_str()), value); }
	void setUniform(const std::string &name, const Vector2i value) const { glUniform2i(glGetUniformLocation(_program, name.c_str()), value[0], value[1]); }
	void setUniform(const std::string &name, const Vector3i value) const { glUniform3i(glGetUniformLocation(_program, name.c_str()), value[0], value[1], value[2]); }
	void setUniform(const std::string &name, const Vector4i value) const { glUniform4i(glGetUniformLocation(_program, name.c_str()), value[0], value[1], value[2], value[3]); }

protected:

	void checkCompileErrors(GLuint object, const std::string &type) const;
};

}
