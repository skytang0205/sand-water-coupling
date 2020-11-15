#pragma once

#include "Utilities/Types.h"

#include <fmt/core.h>

#include <iostream>

namespace PhysX {

class SparseSolver
{
public:

	SparseSolver() = default;
	virtual ~SparseSolver() = default;

	virtual void solve(
		const SparseMatrixr &A,
		Eigen::Ref<VectorXr, Eigen::Aligned> x,
		const Eigen::Ref<const VectorXr, Eigen::Aligned> &b,
		const int max_iterations = 999,
		const real tolerance = real(1e-8)) const = 0;
};

class CgSolver final : public SparseSolver
{
public:

	CgSolver() = default;
	virtual ~CgSolver() = default;

	virtual void solve(
		const SparseMatrixr &A,
		Eigen::Ref<VectorXr, Eigen::Aligned> x,
		const Eigen::Ref<const VectorXr, Eigen::Aligned> &b,
		const int max_iterations = 999,
		const real tolerance = real(1e-8)) const override;
};

class IcPCgSolver final : public SparseSolver
{
public:

	IcPCgSolver() = default;
	virtual ~IcPCgSolver() = default;

	virtual void solve(
		const SparseMatrixr &A,
		Eigen::Ref<VectorXr, Eigen::Aligned> x,
		const Eigen::Ref<const VectorXr, Eigen::Aligned> &b,
		const int max_iterations = 999,
		const real tolerance = real(1e-8)) const override;
};

}
