#include "MaterialPointSubstances.h"

#include "Geometries/ImplicitSurface.h"
#include "Utilities/Constants.h"

#include <fmt/core.h>

#include <numbers>

namespace PhysX {

template <int Dim>
Matrix<Dim, real> MaterialPointSubstances<Dim>::Substance::computeCauchyStressTensor(int idx)
{
	real lambda = _lameLambda;
	real mu = _lameMu;
	if (plastic()) {
		const real ratio = std::exp(_hardeningCoeff * (1 - plasticJacobians[idx]));
		lambda *= ratio;
		mu *= ratio;
	}

	MatrixDr &matF = deformationGradients[idx];

	Eigen::JacobiSVD<MatrixDr> svd(matF, Eigen::ComputeFullU | Eigen::ComputeFullV);
	VectorDr svdS = svd.singularValues();
	const MatrixDr svdU = svd.matrixU();
	const MatrixDr svdV = svd.matrixV();

	if (plastic()) {
		// Handle plasticity.
		plasticJacobians[idx] *= svdS.prod();
		svdS = svdS.cwiseMax(_plasticLowerBound).cwiseMin(_plasticUpperBound);
	}
	const real jacobian = svdS.prod();

	if (plastic()) {
		plasticJacobians[idx] /= jacobian;
		matF = svdU * svdS.asDiagonal() * svdV.transpose();
	}

	if (!mu) {
		// For fluid, reset deformation gradient to avoid numerical instability.
		if constexpr (Dim == 2) matF = MatrixDr::Identity() * std::sqrt(jacobian);
		else matF = MatrixDr::Identity() * std::cbrt(jacobian);
	}

	// compute by the fixed corotated constitutive model.
	return 2 * mu * (matF - svdU * svdV.transpose()) * matF.transpose() + MatrixDr::Identity() * lambda * jacobian * (jacobian - 1);
}

template <int Dim>
void MaterialPointSubstances<Dim>::Substance::reinitialize()
{
	velocities.resize(&particles);
	velocities.setZero();

	velocityDerivatives.resize(&particles);
	velocityDerivatives.setZero();

	deformationGradients.resize(&particles);
	deformationGradients.setConstant(MatrixDr::Identity());

	if (plastic()) {
		plasticJacobians.resize(&particles);
		plasticJacobians.setConstant(1);
	}
}

template <int Dim>
MaterialPointSubstances<Dim>::MaterialPointSubstances(const StaggeredGrid<Dim> &grid) :
	_grid(grid)
{
	_colliders.push_back(
		std::make_unique<StaticCollider<Dim>>(
			std::make_unique<ComplementarySurface<Dim>>(
				std::make_unique<ImplicitBox<Dim>>(
					_grid.domainOrigin() + VectorDr::Ones() * _grid.spacing() / 2,
					_grid.domainLengths() - VectorDr::Ones() * _grid.spacing()))));
}

template <int Dim>
void MaterialPointSubstances<Dim>::writeDescription(YAML::Node &root) const
{
	for (const auto &substance : _substances) {
		YAML::Node node;
		node["name"] = substance.name();
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["material"]["diffuse_albedo"] = substance.color();
		node["indexed"] = false;
		root["objects"].push_back(node);
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	for (const auto &substance : _substances) {
		std::ofstream fout(fmt::format("{}/{}.mesh", frameDir, substance.name()), std::ios::binary);
		IO::writeValue(fout, uint(substance.particles.size()));
		substance.particles.forEach([&](const int i) {
			IO::writeValue(fout, substance.particles.positions[i].template cast<float>().eval());
		});
		if constexpr (Dim == 3) {
			substance.particles.forEach([&](const int i) {
				IO::writeValue(fout, VectorDr::Unit(2).eval());
			});
		}
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::saveFrame(const std::string &frameDir) const
{
	{ // Save velocity.
		std::ofstream fout(frameDir + "/velocity.sav", std::ios::binary);
		_velocity.save(fout);
	}
	// Save substances.
	for (const auto &substance : _substances) {
		std::ofstream fout(fmt::format("{}/{}.sav", frameDir, substance.name()), std::ios::binary);
		IO::writeValue(fout, uint(substance.particles.size()));
		substance.particles.positions.save(fout);
		substance.deformationGradients.save(fout);
		if (substance.plastic())
			substance.plasticJacobians.save(fout);
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::loadFrame(const std::string &frameDir)
{
	{ // Load velocity.
		std::ifstream fin(frameDir + "/velocity.sav", std::ios::binary);
		_velocity.load(fin);
	}
	// Load substances.
	for (auto &substance : _substances) {
		std::ifstream fin(fmt::format("{}/{}.sav", frameDir, substance.name()), std::ios::binary);
		uint particlesCnt;
		IO::readValue(fin, particlesCnt);
		substance.particles.resize(particlesCnt);
		substance.particles.positions.load(fin);

		substance.reinitialize();

		substance.deformationGradients.load(fin);
		if (substance.plastic())
			substance.plasticJacobians.load(fin);
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::initialize()
{
	for (auto &substance : _substances)
		substance.reinitialize();
}

template <int Dim>
void MaterialPointSubstances<Dim>::advance(const real dt)
{
	transferFromGridToParticles(dt);
	applyLagrangianForces(dt);
	transferFromParticlesToGrid(dt);
	applyEulerianForces(dt);
}

template <int Dim>
void MaterialPointSubstances<Dim>::applyEulerianForces(const real dt)
{
	_velocity.parallelForEach([&](const VectorDi &node) {
		if (!_mass[node]) return;
		if (_enableGravity)
			_velocity[node][1] -= kGravity * dt;

		VectorDr tempPos = _velocity.position(node);
		// Resolve collisions.
		for (const auto &collider : _colliders)
			collider->collide(tempPos, _velocity[node]);
	});
}

template <int Dim>
void MaterialPointSubstances<Dim>::transferFromGridToParticles(const real dt)
{
	for (auto &substance : _substances) {
		// Clear particle attributes.
		substance.velocities.setZero();
		substance.velocityDerivatives.setZero();

		substance.particles.parallelForEach([&](const int i) {
			VectorDr &pos = substance.particles.positions[i];
			VectorDr &vel = substance.velocities[i];
			MatrixDr &velDrv = substance.velocityDerivatives[i];
			MatrixDr &defmGrad = substance.deformationGradients[i];
			// Transfer into velocities and velocity derivatives.
			for (const auto [node, weight] : _velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity.position(node) - pos;
				vel += weight * _velocity[node];
				velDrv += _velocity[node] * deltaPos.transpose() * 4 * _velocity.invSpacing() * weight;
			}
			// Advect particles.
			pos += vel * dt;
			defmGrad = (MatrixDr::Identity() + velDrv * dt) * defmGrad;
		});
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::transferFromParticlesToGrid(const real dt)
{
	_velocity.setZero();
	_mass.setZero();

	for (auto &substance : _substances) {
		const real mass = substance.particles.mass();
		const real stressCoeff = dt * 4 * _velocity.invSpacing() * _velocity.invSpacing() * mass / substance.density();

		substance.particles.forEach([&](const int i) {
			const VectorDr pos = substance.particles.positions[i];
			const VectorDr vel = substance.velocities[i];
			const MatrixDr velDrv = substance.velocityDerivatives[i];
			const MatrixDr stress = substance.computeCauchyStressTensor(i) * stressCoeff;
			// Transfer into velocity.
			for (const auto [node, weight] : _velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity.position(node) - pos;
				_velocity[node] += (vel * mass + (velDrv * mass + stress) * deltaPos) * weight;
				_mass[node] += substance.particles.mass() * weight;
			}
		});
	}

	_velocity.parallelForEach([&](const VectorDi &node) {
		if (_mass[node])
			_velocity[node] /= _mass[node];
	});
}

template <int Dim>
void MaterialPointSubstances<Dim>::sampleParticlesInsideSurface(Substance &substance, const Surface<Dim> &surface, const int particlesCntPerSubcell)
{
	substance.particles.clear();

	const int particlesCntPerCell = (1 << Dim) * particlesCntPerSubcell;
	const real dx = _grid.spacing();
	const real radius = dx * real(1.1) / real(std::numbers::sqrt2);
	_grid.forEachCell([&](const VectorDi &cell) {
		const VectorDr centerPos = _grid.cellCenter(cell);
		for (int i = 0; i < particlesCntPerCell; i++) {
			const VectorDr pos = centerPos + VectorDr::Random() * dx / 2;
			if (surface.signedDistance(pos) <= -radius)
				substance.particles.add(pos);
		}
	});

	substance.particles.setMass(substance.density() * std::pow(dx, Dim) / particlesCntPerCell);
}

template class MaterialPointSubstances<2>;
template class MaterialPointSubstances<3>;

}
