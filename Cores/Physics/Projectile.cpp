#include "Projectile.h"

#include "Utilities/Constants.h"
#include "Utilities/IO.h"
#include "Utilities/Yaml.h"

namespace PhysX {

template<int Dim>
void Projectile<Dim>::writeDescription(std::ofstream &fout) const
{
	YAML::Node root;
	root["dimension"] = Dim;
	// Description of ground.
	{
		YAML::Node node;
		node["name"] = "ground";
		node["data_mode"] = "static";
		node["primitive_type"] = "line_list";
		node["indexed"] = false;
		node["enable_color_map"] = false;
		node["diffuse_albedo"] = Vector4f(0, 0, 0, 1);
		root["objects"].push_back(node);
	}
	// Description of ball.
	{
		YAML::Node node;
		node["name"] = "ball";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["indexed"] = false;
		node["enable_color_map"] = false;
		node["diffuse_albedo"] = Vector4f(0, 0, 1, 1);
		root["objects"].push_back(node);
	}
	fout << root << std::endl;
}

template<int Dim>
void Projectile<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	{ // Write ball.
		std::ofstream fout(frameDir + "/ball.mesh", std::ios::binary);
		IO::writeValue(fout, uint(1));
		const VectorDf pos = _position.cast<float>();
		IO::writeValue(fout, pos);
	}
	if (staticDraw) {
		// Write ground.
		std::ofstream fout(frameDir + "/ground.mesh", std::ios::binary);
		IO::writeValue(fout, uint(2));
		IO::writeValue(fout, (VectorDf::Unit(0) * -15.0f).eval());
		IO::writeValue(fout, (VectorDf::Unit(0) * +15.0f).eval());
	}
}

template<int Dim>
void Projectile<Dim>::saveFrame(const std::string &frameDir) const
{
	std::ofstream fout(frameDir + "/state.sav", std::ios::binary);
	IO::writeValue(fout, _position);
	IO::writeValue(fout, _velocity);
}

template<int Dim>
void Projectile<Dim>::loadFrame(const std::string &frameDir)
{
	std::ifstream fin(frameDir + "/state.sav", std::ios::binary);
	IO::readValue(fin, _position);
	IO::readValue(fin, _velocity);
}

template <int Dim>
void Projectile<Dim>::advance(const real dt)
{
	_position += _velocity * dt - VectorDr::Unit(1) * kGravity * dt * dt * real(0.5);
	_velocity.y() -= kGravity * dt;
	_time += dt;
}

template class Projectile<2>;
template class Projectile<3>;

}
