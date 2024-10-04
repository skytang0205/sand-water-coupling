#pragma once

#include "Collider.h"

namespace Pivot {
	class Pressure {
	public:
		explicit Pressure(StaggeredGrid const &sgrid);

	public:
		void Project(SGridData<double> &velocity, GridData<double> const &levelSet, Collider const &collider);
		void Correct(std::vector<Particle> &particles, GridData<double> const &levelSet, Collider const &collider);

		void InitRestDensity(int numPartPerCell, std::vector<Particle> const &particles);

	private:
		void BuildProjectionMatrix(SGridData<double> const &velocity, GridData<double> const &levelSet, Collider const &collider);
		void BuildCorrectionMatrix(std::vector<Particle> const &particles, Collider const &collider);

		void SetUnKnowns(GridData<double> const &levelSet);

		void SolveLinearSystem();

		void ApplyProjection(SGridData<double> &velocity, GridData<double> const &levelSet, Collider const &collider);
		void ApplyCorrection(std::vector<Particle> &particles, Collider const &collider);

	private:
		GridData<int>    m_Grid2Mat;
		std::vector<int> m_Mat2Grid;

		SparseMatrix<double, RowMajor> m_MatL; // A matrix of the Laplacian operator

		VectorXd m_RdP;  // reduced pressure
		VectorXd m_Rhs;

		int               m_NumPartPerCell;
		GridData<double>  m_RestDensity;
		SGridData<double> m_DeltaPos;
	};
}
