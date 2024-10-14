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
		enum class Algorithm { none, alg0, alg1, alg2, alg3};

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
		void MoveDEMParticlesSplit(double ddt, double dt);
		void ApplyDEMForces(double ddt);
		void CalCoupling(double ddt, double dt);
		void TransferCouplingForces(double ddt, double dt);
	
	private:
		// Parameters for basic simulation
		double m_Time;
		Scene  m_Scene;
		// Grid-based data structures
		StaggeredGrid         m_SGrid;
		Collider              m_Collider;
		Pressure              m_Pressure;
		SGridData<double>     m_Velocity;
		SGridData<double> 	  m_ConvectVelocity;
		GridData<double>      m_LevelSet;
		Contour               m_Contour;
		SGridData<double>     m_VelDiff; // Used for FLIP
		// Particle-based data structures
		std::vector<Particle> m_ColliderParticles;
		std::vector<Particle> m_Particles;
		// Parameters for the scheme
		Scheme m_Scheme            		= Scheme::PIC;
		Algorithm m_CouplingAlgorithm   = Algorithm::alg1;
		double m_BlendingFactor    		= 0.95; // Used for FLIP
		double m_LastDeltaTime          = 1000000; //very large
		double m_ViscosityCoeff         = 0.0001;
		// Parameters for particles
		int    m_SeedingSubFactor  		= 3;
		double m_ParticleRadFactor 		= 1.01 * std::numbers::sqrt2 / 2;
		// Boolean configurations
		bool m_DensityCorrectionEnabled = true;
		bool m_SmoothSurfaceEnabled     = true;
		bool m_GravityEnabled           = true;

		//DEM structure

		double m_ParticleRadius;
		double m_SupportRadius;
		double m_ParticleMass;
		double m_ParticleVolume;

		DEMForce m_DEMForce;

		GridData<std::vector<Particle>> m_DEMGrid;

		std::vector<Particle> m_DEMParticles;

		double m_DEMDensity = 1;

		// coupling structure
		SGridData<double>     m_CouplingForce;

	};
}
