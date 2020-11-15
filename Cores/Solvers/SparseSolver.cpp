#include "SparseSolver.h"

namespace PhysX {

void CgSolver::solve(const SparseMatrixr &A, Eigen::Ref<VectorXr, Eigen::Aligned> x, const Eigen::Ref<const VectorXr, Eigen::Aligned> &b, const int max_iterations, const real tolerance) const
{
	Eigen::ConjugateGradient<SparseMatrixr, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> solver(A);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Error: [CgSolver] failed to factorize matrix." << std::endl;
		std::exit(-1);
	}
	solver.setMaxIterations(max_iterations);
	solver.setTolerance(std::max(tolerance, tolerance / b.norm()));
	x = solver.solveWithGuess(b, x);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Error: [CgSolver] failed to solve linear system." << std::endl;
		std::exit(-1);
	}
	std::cout << fmt::format("{:>6} iters", solver.iterations());
}

void IcPCgSolver::solve(const SparseMatrixr &A, Eigen::Ref<VectorXr, Eigen::Aligned> x, const Eigen::Ref<const VectorXr, Eigen::Aligned> &b, const int max_iterations, const real tolerance) const
{
	Eigen::ConjugateGradient<SparseMatrixr, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<real>> solver(A);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Error: [IcPCgSolver] failed to factorize matrix." << std::endl;
		std::exit(-1);
	}
	solver.setMaxIterations(max_iterations);
	solver.setTolerance(std::max(tolerance, tolerance / b.norm()));
	x = solver.solveWithGuess(b, x);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Error: [IcPCgSolver] failed to solve linear system." << std::endl;
		std::exit(-1);
	}
	std::cout << fmt::format("{:>6} iters", solver.iterations());
}

}
