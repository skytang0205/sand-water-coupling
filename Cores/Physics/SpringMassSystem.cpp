#include "SpringMassSystem.h"

#include "Utilities/Constants.h"

namespace PhysX {

template <int Dim>
SpringMassSystem<Dim>::SpringMassSystem() :
	_solver(std::make_unique<IcPCgSolver>())
{ }

template <int Dim>
void SpringMassSystem<Dim>::writeDescription(YAML::Node &root) const
{
	{ // Description of particles.
		YAML::Node node;
		node["name"] = "particles";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["indexed"] = false;
		node["material"]["diffuse_albedo"] = Vector4f(1, 0, 0, 1);
		root["objects"].push_back(node);
	}
	{ // Description of springs.
		YAML::Node node;
		node["name"] = "springs";
		node["data_mode"] = "semi-dynamic";
		node["primitive_type"] = "line_list";
		node["indexed"] = true;
		node["material"]["diffuse_albedo"] = Vector4f(0, 0, 1, 1);
		root["objects"].push_back(node);
	}
}

template <int Dim>
void SpringMassSystem<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	{ // Write particles.
		std::ofstream fout(frameDir + "/particles.mesh", std::ios::binary);
		IO::writeValue(fout, uint(_particles.size()));
		_particles.forEach([&](const int i) {
			IO::writeValue(fout, _particles.positions[i].template cast<float>().eval());
		});
		if constexpr (Dim == 3) {
			_particles.forEach([&](const int i) {
				IO::writeValue(fout, VectorDf::Unit(2).eval());
			});
		}
	}
	{ // Write springs.
		std::ofstream fout(frameDir + "/springs.mesh", std::ios::binary);
		IO::writeValue(fout, uint(_particles.size()));
		_particles.forEach([&](const int i) {
			IO::writeValue(fout, _particles.positions[i].template cast<float>().eval());
		});
		if constexpr (Dim == 3) {
			_particles.forEach([&](const int i) {
				IO::writeValue(fout, VectorDf::Unit(2).eval());
			});
		}
		if (staticDraw) {
			IO::writeValue(fout, 2 * uint(_springs.size()));
			for (const auto &spring : _springs) {
				IO::writeValue(fout, uint(spring.pid0));
				IO::writeValue(fout, uint(spring.pid1));
			}
		}
	}
}

template <int Dim>
void SpringMassSystem<Dim>::saveFrame(const std::string &frameDir) const
{
	{ // Save particles.
		std::ofstream fout(frameDir + "/particles.sav", std::ios::binary);
		_particles.positions.save(fout);
		_velocities.save(fout);
	}
}

template <int Dim>
void SpringMassSystem<Dim>::loadFrame(const std::string &frameDir)
{
	{ // Load particles.
		std::ifstream fin(frameDir + "/particles.sav", std::ios::binary);
		_particles.positions.load(fin);
		_velocities.load(fin);
	}
	reinitializeParticlesBasedData();
}

template <int Dim>
void SpringMassSystem<Dim>::initialize()
{
	reinitializeParticlesBasedData();
}

template <int Dim>
void SpringMassSystem<Dim>::advance(const real dt)
{
	// Advance by the backward Euler method.

	calculateAccelerations();
	buildAndSolveLinearSystem(dt);

	moveParticles(dt);
}

template <int Dim>
void SpringMassSystem<Dim>::calculateAccelerations()
{
	_accelerations.setZero();

	applyElasticForce();
	applyExternalForces();

	// Constrain degrees of freedom.
	for (const int dof : _constrainedDofs) {
		*(reinterpret_cast<real *>(_accelerations.data()) + dof) = 0;
	}
}

template <int Dim>
void SpringMassSystem<Dim>::reinitializeParticlesBasedData()
{
	_accelerations.resize(&_particles);
}

template <int Dim>
void SpringMassSystem<Dim>::buildAndSolveLinearSystem(const real dt)
{
	_matBackwardEuler.resize(_particles.size() * Dim, _particles.size() * Dim);
	_rhsBackwardEuler.resize(_particles.size() * Dim);
	_coeffBackwardEuler.clear();

	const auto addSparseBlockValue = [&](const int pid0, const int pid1, const MatrixDr &block) {
		const int offset0 = pid0 * Dim;
		const int offset1 = pid1 * Dim;
		for (int i = 0; i < Dim; i++) {
			if (_constrainedDofs.contains(offset0 + i)) continue;
			for (int j = 0; j < Dim; j++) {
				if (_constrainedDofs.contains(offset1 + j)) continue;
				_coeffBackwardEuler.push_back(Tripletr(offset0 + i, offset1 + j, block(i, j) * _particles.invMass()));
			}
		}
	};
	const auto asVectorXr = [&](ParticlesBasedVectorData<Dim> &vectorData) {
		return Eigen::Map<VectorXr, Eigen::Aligned>(reinterpret_cast<real *>(vectorData.data()), vectorData.size() * Dim);
	};

	// Set diagonal blocks.
	for (int pid = 0; pid < _particles.size(); pid++)
		addSparseBlockValue(pid, pid, MatrixDr::Identity() * _particles.mass());

	// Set Jacobian of accelaration to velocity.
	for (const auto &spring : _springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = _particles.positions[pid1] - _particles.positions[pid0];
		const VectorDr e01 = r01.normalized();
		const MatrixDr block = spring.dampingCoeff * e01 * e01.transpose() * dt;
		addSparseBlockValue(pid0, pid0, -block);
		addSparseBlockValue(pid0, pid1, block);
		addSparseBlockValue(pid1, pid0, block);
		addSparseBlockValue(pid1, pid1, -block);
	}

	_matBackwardEuler.setFromTriplets(_coeffBackwardEuler.begin(), _coeffBackwardEuler.end());
	_rhsBackwardEuler = _matBackwardEuler * asVectorXr(_velocities) + asVectorXr(_accelerations) * dt;

	// Set Jacobian of acceleration to position.
	for (const auto &spring : _springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = _particles.positions[pid1] - _particles.positions[pid0];
		const double length = r01.norm();
		const VectorDr e01 = r01.normalized();
		const MatrixDr block = spring.stiffnessCoeff * ((spring.restLength / length - 1) * MatrixDr::Identity() - spring.restLength / length * e01 * e01.transpose()) * dt * dt;
		addSparseBlockValue(pid0, pid0, -block);
		addSparseBlockValue(pid0, pid1, block);
		addSparseBlockValue(pid1, pid0, block);
		addSparseBlockValue(pid1, pid1, -block);
	}

	_matBackwardEuler.setFromTriplets(_coeffBackwardEuler.begin(), _coeffBackwardEuler.end());
	_solver->solve(_matBackwardEuler, asVectorXr(_velocities), _rhsBackwardEuler);
}

template <int Dim>
void SpringMassSystem<Dim>::moveParticles(const real dt)
{
	_particles.parallelForEach([&](const int i) {
		_particles.positions[i] += _velocities[i] * dt;
	});

	// Resolve collisions.
	for (const auto &collider : _colliders) {
		collider->collide(_particles, _velocities);
	}
}

template <int Dim>
void SpringMassSystem<Dim>::applyElasticForce()
{
	for (const auto &spring : _springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = _particles.positions[pid1] - _particles.positions[pid0];
		const VectorDr v01 = _velocities[pid1] - _velocities[pid0];
		const real length = r01.norm();
		const VectorDr e01 = r01.normalized();
		const VectorDr force = e01 * ((length - spring.restLength) * spring.stiffnessCoeff + e01.dot(v01) * spring.dampingCoeff);
		_accelerations[pid0] += force * _particles.invMass();
		_accelerations[pid1] -= force * _particles.invMass();
	}
}

template <int Dim>
void SpringMassSystem<Dim>::applyExternalForces()
{
	if (_enableGravity) {
		_particles.parallelForEach([&](const int i) {
			_accelerations[i][1] -= kGravity;
		});
	}
}

template class SpringMassSystem<2>;
template class SpringMassSystem<3>;

}
