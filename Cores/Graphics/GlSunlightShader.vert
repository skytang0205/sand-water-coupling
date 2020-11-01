R"(
#version 450 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

out vec3 vertPos;
out vec3 vertNormal;

uniform mat4 uWorld;

layout (std140, binding = 0) uniform Matrices
{
	mat4 uProjView;
};

void main()
{
	vec4 pos = uWorld * vec4(aPos, 1.0);
	vertPos = pos.xyz;
	vertNormal = mat3(uWorld) * aNormal;
	gl_Position = uProjView * vec4(vertPos, 1.0);
}
)"
