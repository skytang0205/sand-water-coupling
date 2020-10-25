const char *textVsCode = R"
#version 330 core
layout (location = 0) in vec4 aVert;
out vec2 texCoords;

uniform mat4 uProj;

void main()
{
	gl_Position = uProj * vec4(aVert.xy, 0.0, 1.0);
	texCoords = aVert.zw;
}
";

const char *textFsCode = R"
#version 330 core
in vec2 texCoords;
out vec4 fragColor;

uniform sampler2D uText;
uniform vec3 uTextColor;

void main()
{    
    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(uText, texCoords).r);
    fragColor = vec4(uTextColor, 1.0) * sampled;
}
";
