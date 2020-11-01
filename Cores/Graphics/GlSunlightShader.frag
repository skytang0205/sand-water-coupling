R"(
#version 450 core

in vec3 vertPos;
in vec3 vertNormal;

out vec4 fragColor;

void main()
{
	vec3 color = uDiffuseAlbedo.rgb;
	vec3 normal = normalize(vertNormal);
	// Ambient.
	vec3 ambient = uAmbientLight * color;
	// Diffuse.
	vec3 lightDir = -uSunlightDir;
	vec3 diffuse = max(dot(lightDir, normal), 0.0) * color;
	// Specular.
	vec3 viewDir = normalize(uViewPos - vertPos);
	vec3 reflectDir = reflect(-lightDir, normal);
	vec3 halfwayDir = normalize(lightDir + viewDir);
	vec3 specular = pow(max(dot(normal, halfwayDir), 0.0), shininess * 256.0) * FresnelTerm?;

	fragColor = vec4(ambient + diffuse + specular, uDiffuseAlbedo.a);
}
)"
