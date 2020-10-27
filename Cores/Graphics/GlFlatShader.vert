R"(
#version 450 core
layout (location = 0) in vec3 aPos;

uniform mat4 uWorld;

void main()
{
    gl_Position = uWorld * vec4(aPos, 1.0);
}
)"
