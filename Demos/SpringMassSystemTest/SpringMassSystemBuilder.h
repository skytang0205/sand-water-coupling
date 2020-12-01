#pragma once

#include "Physics/SpringMassSystem.h"

namespace PhysX {

class SpringMassSystemBuilder final
{
public:

	template <int Dim>
	static std::unique_ptr<SpringMassSystem<Dim>> build(const int option)
	{
		switch (option) {
		case 0:
			return buildCase0<Dim>();
		case 1:
			return buildCase1<Dim>();
		default:
			reportError("invalid option");
			return nullptr;
		}
	}

protected:

	template <int Dim>
	static std::unique_ptr<SpringMassSystem<Dim>> buildCase0()
	{
		DECLARE_DIM_TYPES(Dim)
		auto smSystem = std::make_unique<SpringMassSystem<Dim>>();
		smSystem->_enableGravity = true;

		smSystem->_masses.add(VectorDr::Ones() * std::numeric_limits<real>::infinity());
		smSystem->_positions.add(VectorDr::Zero());
		smSystem->_velocities.add(VectorDr::Zero());

		smSystem->_masses.add(VectorDr::Ones());
		smSystem->_positions.add(VectorDr::Unit(0));
		smSystem->_velocities.add(VectorDr::Zero());

		smSystem->_springs.push_back({ 0, 1, real(.5), 100, 0 });

		return smSystem;
	}

	template <int Dim>
	static std::unique_ptr<SpringMassSystem<Dim>> buildCase1()
	{
		DECLARE_DIM_TYPES(Dim)
		auto smSystem = std::make_unique<SpringMassSystem<Dim>>();
		for (int i = 0; i < 10; i++) {
			if (i)
				smSystem->_masses.add(VectorDr::Ones());
			else
				smSystem->_masses.add(VectorDr::Ones() * std::numeric_limits<real>::infinity());
			smSystem->_positions.add(VectorDr::Unit(0) * i / 4);
			smSystem->_velocities.add(VectorDr::Zero());
		}
		for (int i = 0; i < 9; i++)
			smSystem->_springs.push_back({ i, i + 1, real(.25), 1000, 0 });
		return smSystem;
	}

	static void reportError(const std::string &msg);
};

}
