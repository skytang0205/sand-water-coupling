#include "MaterialPointSubstances.h"

#include "Geometries/ImplicitSurface.h"
#include "Utilities/Constants.h"

#include <fmt/core.h>

#include <numbers>

namespace PhysX {

template <int Dim>
MaterialPointSubstances<Dim>::MaterialPointSubstances(const StaggeredGrid<Dim> &grid) :
	_grid(grid),
	_velocity(_grid.nodeGrid()),
	_mass(_grid.nodeGrid()),
	_domainBoundary(
		std::make_unique<ComplementarySurface<Dim>>(
			std::make_unique<ImplicitBox<Dim>>(
				_grid.domainOrigin(), _grid.domainLengths()))),
	_integrator(std::make_unique<MpSymplecticEulerIntegrator<Dim>>())
{ }

template <int Dim>
void MaterialPointSubstances<Dim>::writeDescription(YAML::Node &root) const
{
	if constexpr (Dim == 2) { // Description of velocity.
		YAML::Node node;
		node["name"] = "velocity";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "line_list";
		node["indexed"] = false;
		node["color_map"]["enabled"] = true;
		root["objects"].push_back(node);
	}
	for (const auto &substance : _substances) {
		YAML::Node node;
		node["name"] = substance->name();
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["indexed"] = false;
		node["material"]["diffuse_albedo"] = substance->color();
		root["objects"].push_back(node);
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	if constexpr (Dim == 2) { // Write velocity.
		std::ofstream fout(frameDir + "/velocity.mesh", std::ios::binary);
		IO::writeValue(fout, uint(2 * _grid.nodeCount()));
		_grid.forEachNode([&](const VectorDi &node) {
			const VectorDr pos = _grid.nodeCenter(node);
			const VectorDr dir = _velocity[node].normalized() * _grid.spacing() * std::sqrt(real(Dim)) / 2;
			IO::writeValue(fout, pos.template cast<float>().eval());
			IO::writeValue(fout, (pos + dir).template cast<float>().eval());
		});
		_grid.forEachNode([&](const VectorDi &node) {
			const float vel = float(_velocity[node].norm());
			IO::writeValue(fout, vel);
			IO::writeValue(fout, vel);
		});
	}
	for (const auto &substance : _substances) {
		std::ofstream fout(fmt::format("{}/{}.mesh", frameDir, substance->name()), std::ios::binary);
		substance->write(fout);
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
		std::ofstream fout(fmt::format("{}/{}.sav", frameDir, substance->name()), std::ios::binary);
		substance->save(fout);
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
		std::ifstream fin(fmt::format("{}/{}.sav", frameDir, substance->name()), std::ios::binary);
		substance->load(fin);
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::initialize()
{
	for (auto &substance : _substances)
		substance->reinitialize();
}

template <int Dim>
void MaterialPointSubstances<Dim>::advance(const real dt)
{
	transferFromGridToParticles(dt);
	moveParticles(dt);
	applyLagrangianForces(dt);
	transferFromParticlesToGrid(dt);
	applyEulerianForces(dt);
	applyElasticForce(dt);
}

template <int Dim>
void MaterialPointSubstances<Dim>::moveParticles(const real dt)
{
	for (auto &substance : _substances) {
		substance->update(dt);
		// Resolve collisions.
		_domainBoundary.collide(substance->particles.positions);
		for (const auto &collider : _colliders)
			collider->collide(substance->particles.positions);
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::applyEulerianForces(const real dt)
{
	_velocity.parallelForEach([&](const VectorDi &node) {
		if (!_mass[node]) return;
		const VectorDr pos = _velocity.position(node);
		VectorDr &vel = _velocity[node];

		if (_enableGravity)
			vel[1] -= kGravity * dt;

		// Resolve collisions.
		if (_grid.isBoundaryNode(node))
			_domainBoundary.resolve(pos, vel);
		for (const auto &collider : _colliders)
			if (collider->detect(pos)) collider->resolve(pos, vel);
	});
}

template <int Dim>
void MaterialPointSubstances<Dim>::applyElasticForce(const real dt)
{
	_integrator->integrate(_velocity, _mass, _substances, dt);
}

template <int Dim>
void MaterialPointSubstances<Dim>::transferFromGridToParticles(const real dt)
{
	const real gradCoeff = 4 * _velocity.invSpacing() * _velocity.invSpacing();
	for (auto &substance : _substances) {
		// Clear particle attributes.
		substance->velocities.setZero();
		substance->velocityDerivatives.setZero();

		substance->particles.parallelForEach([&](const int i) {
			VectorDr &pos = substance->particles.positions[i];
			VectorDr &vel = substance->velocities[i];
			MatrixDr &velDrv = substance->velocityDerivatives[i];
			// Transfer into velocities and velocity derivatives.
			for (const auto [node, weight] : _velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity.position(node) - pos;
				vel += _velocity[node] * weight;
				velDrv += _velocity[node] * deltaPos.transpose() * gradCoeff * weight;
			}
		});
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::transferFromParticlesToGrid(const real dt)
{
	ParticlesBasedData<Dim, MatrixDr> stresses;
	_velocity.setZero();
	_mass.setZero();

	for (auto &substance : _substances) {
		const real mass = substance->particles.mass();
		const real stressCoeff = -dt * 4 * _velocity.invSpacing() * _velocity.invSpacing() * substance->particles.mass() / substance->density();
		substance->computeStressTensors(stresses);

		substance->particles.forEach([&](const int i) {
			const VectorDr pos = substance->particles.positions[i];
			const VectorDr vel = substance->velocities[i];
			const MatrixDr velDrv = substance->velocityDerivatives[i];
			const MatrixDr stress = stresses[i] * stressCoeff;
			// Transfer into velocity and mass.
			for (const auto [node, weight] : _velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity.position(node) - pos;
				_velocity[node] += (vel * mass + (velDrv * mass + stress) * deltaPos) * weight;
				_mass[node] += mass * weight;
			}
		});
	}

	_velocity.parallelForEach([&](const VectorDi &node) {
		if (_mass[node]) _velocity[node] /= _mass[node];
	});
}

template <int Dim>
void MaterialPointSubstances<Dim>::sampleParticlesInsideSurface(MaterialPointSubstance<Dim> *const substance, const Surface<Dim> &surface, const int particlesCntPerSubcell)
{
	substance->particles.clear();

	const int particlesCntPerCell = (1 << Dim) * particlesCntPerSubcell;
	const real dx = _grid.spacing();
	const real radius = dx * real(1.1) / real(std::numbers::sqrt2);
	_grid.forEachCell([&](const VectorDi &cell) {
		const VectorDr centerPos = _grid.cellCenter(cell);
		for (int i = 0; i < particlesCntPerCell; i++) {
			const VectorDr pos = centerPos + VectorDr::Random() * dx / 2;
			if (Surface<Dim>::isInside(surface.signedDistance(pos) + radius))
				substance->particles.add(pos);
		}
	});

	substance->particles.setMass(substance->density() * std::pow(dx, Dim) / particlesCntPerCell);
}

template class MaterialPointSubstances<2>;
template class MaterialPointSubstances<3>;

}
