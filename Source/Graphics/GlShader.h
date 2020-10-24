#pragma once

#include <glad/glad.h>

#include <string>

namespace PhysX {

class GlShader
{
protected:

	GLuint id;

public:

	GlShader(const std::string &vsFileName, const std::string &fsFileName);

	void Use() const { glUseProgram(id); }

	template <typename Type> void SetUniform(const std::string &name, Type value) const;

#define SPECIALIZE_SET_UNIFORM(type, t) \
template <> void SetUniform(const std::string &name, type value) const { glUniform##t(glGetUniformLocation(id, name.c_str()), value); }

	SPECIALIZE_SET_UNIFORM(GLint, 1i)
	SPECIALIZE_SET_UNIFORM(GLfloat, 1f)

#undef SPECIALIZE_SET_UNIFORM

protected:

	void checkCompileErrors(GLuint object, const std::string &type) const;
};

}
