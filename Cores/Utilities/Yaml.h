#pragma once
#pragma warning (disable : 4996)

#include "Utilities/Types.h"

#include <yaml-cpp/yaml.h>

namespace YAML {

template <class Scalar, int Dim>
struct convert<PhysX::Vector<Scalar, Dim>>
{
	static Node encode(const PhysX::Vector<Scalar, Dim> &rhs)
	{
		Node node;
		for (int i = 0; i < rhs.size(); i++) {
			node.push_back(rhs(i));
		}
		node.SetStyle(EmitterStyle::Flow);
		return node;
	}
	static bool decode(const Node &node, PhysX::Vector<Scalar, Dim> &rhs)
	{
		if (!node.IsSequence() || node.size() != rhs.size()) {
			return false;
		}
		for (int i = 0; i < rhs.size(); i++) {
			rhs(i) = node[i].as<Scalar>();
		}
		return true;
	}
};

}
