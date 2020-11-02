R"(
#version 450 core

in vec3 vertPos;
in vec3 vertNormal;

out vec4 fragColor;

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
	float uTotalTime;		// 128
	float uDeltaTime;		// 132
};

vec3 calcSchlickFresnel(const vec3 R0, const vec3 normal, const vec3 lightDir)
{
	float cosIncidentAngle = clamp(dot(normal, lightDir), 0.0, 1.0);
	float f0 = 1.0 - cosIncidentAngle;
	return R0 + (1.0 - R0) * (f0 * f0 * f0 * f0 * f0);
}

vec3 calcDiffuseAndSpecular(const vec3 lightStrength, const vec3 lightDir, const vec3 normal, const vec3 viewDir)
{
	vec3 halfDir = normalize(lightDir + viewDir);
	float m = (1.0 - uRoughness) * 256.0;
	float roughnessFactor = (m + 8.0) * pow(max(dot(halfDir, normal), 0.0), m) / 8.0;
	vec3 fresnelFactor = calcSchlickFresnel(uFresnelR0, halfDir, lightDir);
	vec3 specularAlbedo = fresnelFactor * roughnessFactor;
	// Apply LDR.
	specularAlbedo = specularAlbedo / (specularAlbedo + 1.0);
	return (uDiffuseAlbedo.rgb + specularAlbedo) * lightStrength;
}

vec3 ComputeDirectionalLight(const vec3 normal, const vec3 viewDir)
{
	vec3 lightStrength = uLightStrength * max(dot(uLightDir, normal), 0.0);
	return calcDiffuseAndSpecular(lightStrength, uLightDir, normal, viewDir);
}

void main()
{
	vec3 ambient = uAmbientStrength * uDiffuseAlbedo.rgb;
	vec3 diffspec = ComputeDirectionalLight(normalize(vertNormal), normalize(uViewPos - vertPos));
	fragColor = vec4(ambient + diffspec, uDiffuseAlbedo.a);
}
)"
