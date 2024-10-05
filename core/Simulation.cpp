#include "Simulation.h"

#include "Advection.h"
#include "TriBSpline.h"
#include "TriLerp.h"
#include "CSG.h"
#include "Extrapolation.h"
#include "FiniteDiff.h"
#include "Redistancing.h"

namespace Pivot {
	Simulation::Simulation(StaggeredGrid const &sgrid, double particleRadius) :
		m_SGrid { sgrid },
		m_Collider(m_SGrid),
		m_Pressure(m_SGrid),
		m_Velocity(m_SGrid.GetFaceGrids()),
		m_LevelSet(m_SGrid.GetCellGrid(), std::numeric_limits<double>::infinity()),
		m_Contour(m_SGrid.GetCellGrid()),
		m_VelDiff(m_SGrid.GetFaceGrids()),
		m_DEMGrid(m_SGrid.GetCellGrid()),
		m_DEMForce(particleRadius),
		m_ParticleRadius(particleRadius),
		m_SupportRadius(2 * particleRadius) {
	}

	void Simulation::Describe(YAML::Node &root) const {
		root["Dimension"] = 3;
		root["Radius"] = m_SGrid.GetDomainRadius();
		{ // Description of particles
			YAML::Node node;
			node["Name"] = "particles";
			node["Animated"] = true;
			node["Primitive"] = "Points";
			node["Material"]["Albedo"] = Vector4f(0, 0, 1, 1);
			root["Objects"].push_back(node);
		}
		{ // Description of particles
			YAML::Node node;
			node["Name"] = "DEMparticles";
			node["Animated"] = true;
			node["Primitive"] = "Points";
			node["Material"]["Albedo"] = Vector4f(0, 0, 1, 1);
			root["Objects"].push_back(node);
		}
		{ // Description of collider
			YAML::Node node;
			node["Name"] = "collider";
			node["Primitive"] = "Points";
			node["Material"]["Albedo"] = Vector4f(.5f, .5f, .5f, 1);
			root["Objects"].push_back(node);
		}
		// { // Description of contour
		// 	YAML::Node node;
		// 	node["Name"] = "contour";
		// 	node["Animated"] = true;
		// 	node["Indexed"] = true;
		// 	node["Primitive"] = "Triangles";
		// 	node["Material"]["Albedo"] = Vector4f(0, 0, 1, 1);
		// 	root["Objects"].push_back(node);
		// }
	}

	void Simulation::Export(std::filesystem::path const &dirname, bool initial) const {
		{ // Export particles
			std::ofstream fout(dirname / "particles.out", std::ios::binary);
			IO::Write(fout, static_cast<std::uint32_t>(m_Particles.size()));
			for (auto const &particle : m_Particles) {
				IO::Write(fout, particle.Position.cast<float>().eval());
			}
			for (auto const &particle : m_Particles) {
			}
		}
		{ // Export regular particles
			std::ofstream fout(dirname / "DEMparticles.out", std::ios::binary);
			IO::Write(fout, static_cast<std::uint32_t>(m_DEMParticles.size()));
			for (auto const &p : m_DEMParticles) {
				IO::Write(fout, p.Position.cast<float>().eval());
			}
			// for (auto const &p : m_DEMParticles) {
			// 	IO::Write(fout, static_cast<float>(p.Velocity.norm()));
			// }
			for (auto const &p : m_DEMParticles) {
				IO::Write(fout, static_cast<float>(m_ParticleRadius));
			}
		}
		// { // Export the contour
		// 	std::ofstream fout(dirname / "contour.out", std::ios::binary);
		// 	m_Contour.GetMesh().Export(fout);
		// }
		if (initial) { // Export the collider
			std::ofstream fout(dirname / "collider.out", std::ios::binary);
			std::vector<Vector3f> positions;
			std::vector<Vector3f> normals;
			for (auto const &particle : m_ColliderParticles) {
				if (m_Collider.GetDomainBox().SignedDistanceTo(particle.Position) < -.01 * m_SGrid.GetSpacing()) {
					positions.push_back(particle.Position.cast<float>());
					normals.push_back(TriLerp::Interpolate(m_Collider.GetNormal(), particle.Position).normalized().cast<float>());
				}
			}
			IO::Write(fout, static_cast<std::uint32_t>(positions.size()));
			IO::Write(fout, positions);
			IO::Write(fout, normals);
		}
	}

	void Simulation::Save(std::ostream &out) const {
	}

	void Simulation::Load(std::istream &in) {
	}

	double Simulation::GetCourantTimeStep() const { 
		double maxVel = 0;
		for (auto const &particle : m_DEMParticles) {
			maxVel = std::max(maxVel, particle.Velocity.norm());
		}
		//printf("%lf\n",maxVel);
		maxVel += std::sqrt(9.8 * m_ParticleRadius * 2);
		

		return  std::min(m_SGrid.GetSpacing() / (m_Velocity.GetMaxAbsComponent()), m_ParticleRadius * 2. / maxVel); 
	}

	void Simulation::Initialize() {
		m_Collider.Finish(m_SGrid);
		CSG::Intersect(m_LevelSet, m_Collider.GetDomainBox());
		
		SeedParticles();
		ReconstructLevelSet();
		GenerateContour();

		if (m_DensityCorrectionEnabled)
			m_Pressure.InitRestDensity(m_SeedingSubFactor * m_SeedingSubFactor * m_SeedingSubFactor, m_ColliderParticles);

		m_ParticleMass = m_TargetDensity * m_SupportRadius * m_SupportRadius * m_SupportRadius * std::numbers::sqrt2 * .5;
		CacheNeighborHoods();
	}

	void Simulation::Advance(double deltaTime) {
		
		TransferFromGridToParticles();
		AdvectFields(deltaTime);
		TransferFromParticlesToGrid();

		ApplyBodyForces(deltaTime);
		ApplyDEMForces(deltaTime);
		MoveDEMParticles(deltaTime);

		ProjectVelocity(deltaTime);
		CacheNeighborHoods();
	}

	void Simulation::TransferFromGridToParticles() {
		switch (m_Scheme) {
		case Scheme::PIC:
			tbb::parallel_for_each(m_Particles.begin(), m_Particles.end(), [&](Particle &particle) {
				for (int axis = 0; axis < 3; axis++) {
					particle.Velocity[axis] = TriLerp::Interpolate(m_Velocity[axis], particle.Position);
				}
			});
			break;
		case Scheme::FLIP:
			tbb::parallel_for_each(m_Particles.begin(), m_Particles.end(), [&](Particle &particle) {
				for (int axis = 0; axis < 3; axis++) {
					particle.Velocity[axis] += TriLerp::Interpolate(m_VelDiff[axis], particle.Position);
					particle.Velocity[axis] = m_BlendingFactor * particle.Velocity[axis] + (1 - m_BlendingFactor) * TriLerp::Interpolate(m_Velocity[axis], particle.Position);
				}
			});
			break;
		case Scheme::APIC:
			tbb::parallel_for_each(m_Particles.begin(), m_Particles.end(), [&](Particle &particle) {
				for (int axis = 0; axis < 3; axis++) {
					particle.Velocity[axis] = TriLerp::Interpolate(m_Velocity[axis], particle.Position);
					particle.VelocityDrv[axis].setZero();
					for (auto const [face, weight] : TriLerp::GetGradWtPoints(m_Velocity[axis].GetGrid(), particle.Position)) {
						particle.VelocityDrv[axis] += m_Velocity[axis].At(face) * weight;
					}
				}
			});
			break;
		}
	}

	void Simulation::TransferFromParticlesToGrid() {
		m_Velocity.SetZero();
		SGridData<double> weightSum(m_Velocity.GetGrids());

		if (m_Scheme != Scheme::APIC) {
			for (auto const &particle : m_Particles) {
				for (int axis = 0; axis < 3; axis++) {
					for (auto [face, weight] : TriLerp::GetWtPoints(m_Velocity[axis].GetGrid(), particle.Position)) {
						face = m_Velocity[axis].GetGrid().Clamp(face);
						m_Velocity[axis][face] += particle.Velocity[axis] * weight;
						weightSum[axis][face] += weight;
					}
				}
			}
		} else {
			for (auto const &particle : m_Particles) {
				for (int axis = 0; axis < 3; axis++) {
					for (auto [face, weight] : TriLerp::GetWtPoints(m_Velocity[axis].GetGrid(), particle.Position)) {
						Vector3d const deltaPos = m_Velocity[axis].GetGrid().PositionOf(face) - particle.Position;
						face = m_Velocity[axis].GetGrid().Clamp(face);
						m_Velocity[axis][face] += (particle.Velocity[axis] + particle.VelocityDrv[axis].dot(deltaPos)) * weight;
						weightSum[axis][face] += weight;
					}
				}
			}
		}

		ParallelForEach(m_Velocity.GetGrids(), [&](int axis, Vector3i const &face) {
			if (weightSum[axis][face]) {
				m_Velocity[axis][face] /= weightSum[axis][face];
			}
		});

		Extrapolation::Solve(m_Velocity, 0., 3, [&](int axis, Vector3i const &face) {
			return weightSum[axis][face] != 0;
		});

		if (m_Scheme == Scheme::FLIP) {
			m_VelDiff = m_Velocity;
		}
	}

	void Simulation::AdvectFields(double dt) {
		tbb::parallel_for_each(m_Particles.begin(), m_Particles.end(), [&](Particle &particle) {
			particle.Position = AdvectionSL::Trace<2>(particle.Position, m_Velocity, dt);
		});

		if (m_DensityCorrectionEnabled) {
			ReconstructLevelSet();
			m_Pressure.Correct(m_Particles, m_LevelSet, m_Collider);
		}

		m_Collider.Enforce(m_Particles);

		ReconstructLevelSet();
		GenerateContour();
	}

	void Simulation::ApplyBodyForces(double dt) {
		if (m_GravityEnabled) {
			ParallelForEach(m_Velocity[1].GetGrid(), [&](Vector3i const &face) {
				m_Velocity[1][face] -= 9.8 * dt;
			});
			tbb::parallel_for_each(m_DEMParticles.begin(), m_DEMParticles.end(), [&](Particle &p) {
				p.Velocity[1] -= 9.8 * dt;
			});
		}
	}

	void Simulation::ProjectVelocity(double dt) {
		m_Pressure.Project(m_Velocity, m_LevelSet, m_Collider);
		Extrapolation::Solve(m_Velocity, 0., 6, [&](int axis, Vector3i const &face) {
			Vector3i const cell0 = StaggeredGrid::AdjCellOfFace(axis, face, 0);
			Vector3i const cell1 = StaggeredGrid::AdjCellOfFace(axis, face, 1);
			return m_Collider.GetFraction()[axis][face] < 1 && (m_LevelSet[cell0] <= 0 || m_LevelSet[cell1] <= 0);
		});
		m_Collider.Enforce(m_Velocity);

		if (m_Scheme == Scheme::FLIP) {
			ParallelForEach(m_VelDiff.GetGrids(), [&](int axis, Vector3i const &face) {
				m_VelDiff[axis][face] = m_Velocity[axis][face] - m_VelDiff[axis][face];
			});
		}
	}

	void Simulation::SeedParticles() {
		Extrapolation::Solve(m_LevelSet, 1.5 * m_SGrid.GetSpacing(), 1, [&](Vector3i const &cell) {
			return !m_Collider.IsInside(cell);
		});
		RedistancingFM::Solve(m_LevelSet, 5);

		// A stratified sampling
		auto const rgrid = m_SGrid.CreateRefinedCellGrid(m_SeedingSubFactor);
		double const dx  = m_SGrid.GetSpacing();
		double const rad = m_ParticleRadFactor * dx;

		m_Particles.clear();
		ForEach(rgrid, [&](Vector3i const &cell) {
			Vector3d const pos = rgrid.PositionOf(cell) + Vector3d::Random() * rgrid.GetSpacing() / 2;
			if (TriLerp::Interpolate(m_Collider.LevelSet, pos) <= 0) {
				m_ColliderParticles.push_back({
					.Position = pos,
				});
			} else if (TriLerp::Interpolate(m_LevelSet, pos) + rad <= 0) {
				m_Particles.push_back({
					.Position = pos,
				});
			}
		});
	}

	void Simulation::ReconstructLevelSet() {
		double const dx  = m_SGrid.GetSpacing();
		double const rad = m_ParticleRadFactor * dx;
		m_LevelSet.SetConstant(2 * dx - rad);
		for (auto const &particle : m_Particles) {
			for (auto const cell : TriBSpline<3>::GetPoints(m_LevelSet.GetGrid(), particle.Position)) {
				if (!m_LevelSet.GetGrid().IsValid(cell)) break;
				Vector3d const pos = m_LevelSet.GetGrid().PositionOf(cell);
				m_LevelSet[cell] = std::min(m_LevelSet[cell], (particle.Position - pos).norm() - rad);
			}
		}
		RedistancingFM::Solve(m_LevelSet, 5);

		if (m_SmoothSurfaceEnabled) {
			auto const oldLevelSet = m_LevelSet;
			ParallelForEach(m_LevelSet.GetGrid(), [&](Vector3i const &cell) {
				double mean = 0;
				for (int i = 0; i < Grid::GetNumNeighbors(); i++) {
					Vector3i const nbCell = Grid::NeighborOf(cell, i);
					mean += oldLevelSet.At(nbCell);
				}
				mean /= Grid::GetNumNeighbors();
				if (mean < m_LevelSet[cell]) m_LevelSet[cell] = mean;
			});
			RedistancingFM::Solve(m_LevelSet, 5);
		}
	}

	void Simulation::GenerateContour() {
		auto opLevelSet = m_LevelSet;
		CSG::Except(opLevelSet, m_Collider.GetAuxLevelSet());
		m_Contour.Generate(opLevelSet);
		m_Contour.ComputeVertexInfos();
	}

	void Simulation::CacheNeighborHoods() {
		ParallelForEach(m_DEMGrid.GetGrid(), [&](Vector3i const &cell){
			m_DEMGrid[cell].clear();
		});
		//printf("breakpoint 1\n");
		for(Particle p : m_DEMParticles){
			//printf("breakpoint 3\n");
			Vector3i lower = m_DEMGrid.GetGrid().Clamp(m_DEMGrid.GetGrid().CalcLower<1>(p.Position));
			//printf("breakpoint 4\n");
			//std::lock_guard<std::mutex> lock(m_mutex);
			m_DEMGrid[lower].push_back(p);
		}
		//printf("breakpoint 2\n");
	}

	void Simulation::MoveDEMParticles(double dt) {
		tbb::parallel_for_each(m_DEMParticles.begin(), m_DEMParticles.end(), [&](Particle &p) {
			p.Position += p.Velocity * dt;
		});
		m_Collider.Enforce(m_DEMParticles);
	}

	void Simulation::ApplyDEMForces(double dt) {
		// Note: Do not apply external forces to boundary particles
		tbb::parallel_for_each(m_DEMParticles.begin(), m_DEMParticles.end(), [&](Particle &p) {
				p.Velocity += m_DEMForce.getForceSum(m_DEMGrid, p) * dt / m_ParticleMass;
			});		
	}
}
