#pragma once

#include "Materials/MaterialPointSubstance.h"
#include "Solvers/SparseSolver.h"
#include "Structures/GridBasedData.h"

namespace PhysX {

template <int Dim>
class MaterialPointIntegrator
{
	DECLARE_DIM_TYPES(Dim)

public:

	MaterialPointIntegrator() = default;
	MaterialPointIntegrator(const MaterialPointIntegrator &rhs) = delete;
	MaterialPointIntegrator &operator=(const MaterialPointIntegrator &rhs) = delete;
	virtual ~MaterialPointIntegrator() = default;

	virtual void integrate(
		GridBasedVectorData<Dim> &velocity,
		const GridBasedScalarData<Dim> &mass,
		const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances,
		const real dt,
		const GridBasedData<Dim, uchar> &collided) = 0;
};

template <int Dim>
class MpSymplecticEulerIntegrator : public MaterialPointIntegrator<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	MpSymplecticEulerIntegrator() = default;
	MpSymplecticEulerIntegrator(const MpSymplecticEulerIntegrator &rhs) = delete;
	MpSymplecticEulerIntegrator &operator=(const MpSymplecticEulerIntegrator &rhs) = delete;
	virtual ~MpSymplecticEulerIntegrator() = default;

	virtual void integrate(
		GridBasedVectorData<Dim> &velocity,
		const GridBasedScalarData<Dim> &mass,
		const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances,
		const real dt,
		const GridBasedData<Dim, uchar> &collided) override
	{ }
};

template <int Dim>
class MpSemiImplicitIntegrator : public MaterialPointIntegrator<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::vector<Tripletr> _coefficients;
	SparseMatrixr _matLinearized;
	VectorXr _rhsLinearized;

	std::unique_ptr<SparseSolver> _solver;

public:

	MpSemiImplicitIntegrator();

	MpSemiImplicitIntegrator(const MpSemiImplicitIntegrator &rhs) = delete;
	MpSemiImplicitIntegrator &operator=(const MpSemiImplicitIntegrator &rhs) = delete;
	virtual ~MpSemiImplicitIntegrator() = default;

	virtual void integrate(
		GridBasedVectorData<Dim> &velocity,
		const GridBasedScalarData<Dim> &mass,
		const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances,
		const real dt,
		const GridBasedData<Dim, uchar> &collided) override;
};

}
