#pragma once

#include "Collider.h"
#include "Contour.h"
#include "Pressure.h"
#include "DEMForce.h"

namespace Pivot {
	class Simulation {
	private:
		friend class SimBuilder;
	
	public:
		enum class Scheme { PIC, FLIP, APIC };
		enum class Scene { Falling, BigBall, Slope };

	public:
		explicit Simulation(StaggeredGrid const &sgrid, double particleRadius);

		void Describe(YAML::Node &root) const;
		void Export(std::filesystem::path const &dirname, bool initial = false) const;
		void Save(std::ostream &out) const;
		void Load(std::istream &in);

		void Initialize();
		void Advance(double deltaTime);

		void TransferFromGridToParticles();
		void TransferFromParticlesToGrid();
		void AdvectFields(double dt);
		void ApplyBodyForces(double dt);
		void ProjectVelocity(double dt);

		void SeedParticles();
		void ReconstructLevelSet();
		void GenerateContour();

		void SetTime(double time) { m_Time = time; }
		auto GetTime() const { return m_Time; }

		double GetCourantTimeStep() const;

		void CacheNeighborHoods();

		void MoveDEMParticles(double dt);
		void ApplyDEMForces(double dt);
	
	private:
		// Parameters for basic simulation
		double m_Time;
		Scene  m_Scene;
		// Grid-based data structures
		StaggeredGrid         m_SGrid;
		Collider              m_Collider;
		Pressure              m_Pressure;
		SGridData<double>     m_Velocity;
		GridData<double>      m_LevelSet;
		Contour               m_Contour;
		SGridData<double>     m_VelDiff; // Used for FLIP
		// Particle-based data structures
		std::vector<Particle> m_ColliderParticles;
		std::vector<Particle> m_Particles;
		// Parameters for the scheme
		Scheme m_Scheme            = Scheme::APIC;
		double m_BlendingFactor    = 0.95; // Used for FLIP
		// Parameters for particles
		int    m_SeedingSubFactor  = 3;
		double m_ParticleRadFactor = 1.01 * std::numbers::sqrt2 / 2;
		// Boolean configurations
		bool m_DensityCorrectionEnabled = true;
		bool m_SmoothSurfaceEnabled     = true;
		bool m_GravityEnabled           = true;

		//DEM structure

		double m_ParticleRadius;
		double m_SupportRadius;
		double m_ParticleMass;

		DEMForce m_DEMForce;

		GridData<std::vector<Particle>> m_DEMGrid;

		std::vector<Particle> m_DEMParticles;

		double m_TargetDensity = 2700;
	};
}
