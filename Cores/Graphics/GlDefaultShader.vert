R"(
#version 450 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

out vec3 vertPos;
out vec3 vertNormal;

uniform mat4 uWorld;
uniform vec4 uDiffuseAlbedo;
uniform vec3 uFresnelR0;
uniform float uRoughness;

layout (std140, binding = 0) uniform PassConstants
{
	mat4 uProjView;			// 0
	vec3 uViewPos;			// 64
	vec3 uAmbientStrength;	// 80
	vec3 uLightStrength;	// 96
	vec3 uLightDir;			// 112
	float uTotalTime;		// 124
	float uDeltaTime;		// 128
};

void main()
{
	vec4 pos = uWorld * vec4(aPos, 1.0);
	vertPos = pos.xyz;
	vertNormal = mat3(uWorld) * aNormal;
	gl_Position = uProjView * vec4(vertPos, 1.0);
}
)"
