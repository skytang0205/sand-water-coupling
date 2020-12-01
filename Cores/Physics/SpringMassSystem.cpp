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
		IO::writeValue(fout, uint(_positions.size()));
		_positions.forEach([&](const int i) {
			IO::writeValue(fout, _positions[i].cast<float>().eval());
		});
	}
	{ // Write springs.
		std::ofstream fout(frameDir + "/springs.mesh", std::ios::binary);
		IO::writeValue(fout, uint(_positions.size()));
		_positions.forEach([&](const int i) {
			IO::writeValue(fout, _positions[i].cast<float>().eval());
		});
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
	{ // Save positions.
		std::ofstream fout(frameDir + "/positions.sav", std::ios::binary);
		_positions.save(fout);
	}
	{ // Save velocities.
		std::ofstream fout(frameDir + "/velocities.sav", std::ios::binary);
		_velocities.save(fout);
	}
}

template <int Dim>
void SpringMassSystem<Dim>::loadFrame(const std::string &frameDir)
{
	{ // Load positions.
		std::ifstream fin(frameDir + "/positions.sav", std::ios::binary);
		_positions.load(fin);
	}
	{ // Load velocities.
		std::ifstream fin(frameDir + "/velocities.sav", std::ios::binary);
		_velocities.load(fin);
	}
	// Reinitialize attributs and reset the sparse matrix.
	reinitializeAttributes();
	resetSparse();
}

template <int Dim>
void SpringMassSystem<Dim>::initialize()
{
	reinitializeAttributes();
	resetSparse();
}

template <int Dim>
void SpringMassSystem<Dim>::advance(const real dt)
{
	// Advance by the backward Euler method.
	_matBackwardEuler.resize(_positions.size() * Dim, _positions.size() * Dim);
	_matBackwardEuler.setZero();

	calculateAccelarations();

	buildAndSolveLinearSystem(dt);

	// Update positions.
	_positions.parallelForEach([&](const int pid) {
		_positions[pid] += _velocities[pid] * dt;
	});
}

template <int Dim>
void SpringMassSystem<Dim>::calculateAccelarations()
{
	_accelerations.setZero();
	applyElasticForces();
	applyExternalForces();
}

template <int Dim>
void SpringMassSystem<Dim>::reinitializeAttributes()
{
	_invMasses.resize(_masses.size());
	_invSqrtMasses.resize(_masses.size());
	_accelerations.resize(_masses.size());

	_masses.parallelForEach([&](const int pid) {
		_invMasses[pid] = _masses[pid].cwiseInverse();
		_invSqrtMasses[pid] = _invMasses[pid].cwiseSqrt();
	});
}

template <int Dim>
void SpringMassSystem<Dim>::resetSparse()
{
	_matBackwardEuler.resize(_masses.size() * Dim, _masses.size() * Dim);
	_rhsBackwardEuler.resize(_masses.size() * Dim);
	std::vector<Tripletr> elements;

	const auto addSparseTriplets = [&](const int pid0, const int pid1) {
		const int offset0 = pid0 * Dim;
		const int offset1 = pid1 * Dim;
		for (int i = 0; i < Dim; i++)
			for (int j = 0; j < Dim; j++)
				elements.push_back(Tripletr(offset0 + i, offset1 + j));
	};

	for (int pid = 0; pid < _positions.size(); pid++)
		addSparseTriplets(pid, pid);
	for (const auto &spring : _springs) {
		addSparseTriplets(spring.pid0, spring.pid1);
		addSparseTriplets(spring.pid1, spring.pid0);
	}
	_matBackwardEuler.setFromTriplets(elements.begin(), elements.end());
}

template <int Dim>
void SpringMassSystem<Dim>::buildAndSolveLinearSystem(const real dt)
{
	const auto addSparseBlockValue = [&](const int pid0, const int pid1, const MatrixDr &block) {
		const size_t offset0 = size_t(pid0) * Dim;
		const size_t offset1 = size_t(pid1) * Dim;
		const MatrixDr reducedBlock = _invSqrtMasses[pid0].asDiagonal() * block * _invSqrtMasses[pid1].asDiagonal();
		for (int i = 0; i < Dim; i++)
			for (int j = 0; j < Dim; j++)
				_matBackwardEuler.coeffRef(offset0 + i, offset1 + j) += reducedBlock(i, j);
	};
	const auto asVectorXr = [&](ParticlesVectorAttribute<Dim> &attr) {
		return Eigen::Map<VectorXr, Eigen::Aligned>(reinterpret_cast<real *>(attr.data()), attr.size() * Dim);
	};

	// Set diagonal blocks.
	for (int pid = 0; pid < _positions.size(); pid++)
		addSparseBlockValue(pid, pid, MatrixDr::Identity());

	// Set Jacobian of accelaration to velocity.
	for (const auto &spring : _springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = _positions[pid1] - _positions[pid0];
		const VectorDr e01 = r01.normalized();
		const MatrixDr block = spring.dampingCoeff * e01 * e01.transpose() * dt;
		addSparseBlockValue(pid0, pid0, -block);
		addSparseBlockValue(pid0, pid1, block);
		addSparseBlockValue(pid1, pid0, block);
		addSparseBlockValue(pid1, pid1, -block);
	}

	_rhsBackwardEuler = _matBackwardEuler * asVectorXr(_velocities) + asVectorXr(_accelerations) * dt;

	// Set Jacobian of acceleration to position.
	for (const auto &spring : _springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = _positions[pid1] - _positions[pid0];
		const double length = r01.norm();
		const VectorDr e01 = r01.normalized();
		const MatrixDr block = spring.stiffnessCoeff * ((spring.restLength / length - 1) * MatrixDr::Identity() - spring.restLength / length * e01 * e01.transpose()) * dt * dt;
		addSparseBlockValue(pid0, pid0, -block);
		addSparseBlockValue(pid0, pid1, block);
		addSparseBlockValue(pid1, pid0, block);
		addSparseBlockValue(pid1, pid1, -block);
	}
	_solver->solve(_matBackwardEuler, asVectorXr(_velocities), _rhsBackwardEuler);
}

template <int Dim>
void SpringMassSystem<Dim>::applyElasticForces()
{
	for (const auto &spring : _springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = _positions[pid1] - _positions[pid0];
		const VectorDr v01 = _velocities[pid1] - _velocities[pid0];
		const real length = r01.norm();
		const VectorDr e01 = r01.normalized();
		const VectorDr force = e01 * ((length - spring.restLength) * spring.stiffnessCoeff + e01.dot(v01) * spring.dampingCoeff);
		_accelerations[pid0] += force.cwiseProduct(_invMasses[pid0]);
		_accelerations[pid1] -= force.cwiseProduct(_invMasses[pid1]);
	}
}

template <int Dim>
void SpringMassSystem<Dim>::applyExternalForces()
{
	if (_enableGravity) {
		_accelerations.parallelForEach([&](const int pid) {
			if (_invMasses[pid][1]) _accelerations[pid][1] -= kGravity;
		});
	}
}

template class SpringMassSystem<2>;
template class SpringMassSystem<3>;

}
