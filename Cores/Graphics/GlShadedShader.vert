R"(
#version 450 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec4 aColor;

out vec4 vertColor;

layout (std140, binding = 0) uniform Matrices
{
	mat4 uProjView;
};

void main()
{
	gl_Position = uProjView * vec4(aPos, 1.0);
	vertColor = aColor;
}
)"
