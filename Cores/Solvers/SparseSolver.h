#pragma once

#include "Utilities/Types.h"

#include <fmt/core.h>

#include <iostream>

namespace PhysX {

class SparseSolver
{
protected:

	static constexpr int _kDefaultMaxIterations = 999;
	static constexpr real _kDefaultTolerance = real(1e-6);

public:

	SparseSolver() = default;
	virtual ~SparseSolver() = default;

	virtual void solve(
		const SparseMatrixr &A,
		Eigen::Ref<VectorXr, Eigen::Aligned> x,
		const Eigen::Ref<const VectorXr, Eigen::Aligned> &b,
		const int maxIterations = _kDefaultMaxIterations,
		const real tolerance = _kDefaultTolerance) const = 0;
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
		const int maxIterations = _kDefaultMaxIterations,
		const real tolerance = _kDefaultTolerance) const override;
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
		const int maxIterations = _kDefaultMaxIterations,
		const real tolerance = _kDefaultTolerance) const override;
};

}
