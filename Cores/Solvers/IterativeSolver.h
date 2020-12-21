#pragma once

#include "Utilities/Types.h"

#include <fmt/core.h>

#include <iostream>

namespace PhysX::IterativeSolver {

template <typename MatrixType>
using CG = Eigen::ConjugateGradient<MatrixType, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner>;

template <typename MatrixType>
using ICPCG = Eigen::ConjugateGradient<MatrixType, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<real>>;

template <typename MatrixType>
using BiCGSTAB = Eigen::BiCGSTAB<MatrixType, Eigen::IdentityPreconditioner>;

template <typename MatrixType, typename Solver = ICPCG<MatrixType>>
inline void solve(const MatrixType &A, Eigen::Ref<VectorXr, Eigen::Aligned> x, const Eigen::Ref<const VectorXr, Eigen::Aligned> &b, const int maxIterations = -1, const real tolerance = real(1e-6))
{
	Solver solver(A);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Error: [IterativeSolver] failed to factorize matrix." << std::endl;
		std::exit(-1);
	}
	if (maxIterations >= 0) solver.setMaxIterations(maxIterations);
	solver.setTolerance(std::max(tolerance, std::numeric_limits<real>::epsilon() / b.norm()));
	x = solver.solveWithGuess(b, x);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Error: [IterativeSolver] failed to solve linear system." << std::endl;
		std::exit(-1);
	}
	std::cout << fmt::format("{:>6} iters", solver.iterations());
}

}
