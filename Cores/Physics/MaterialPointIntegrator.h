#pragma once

#include "Materials/MaterialPointSubstance.h"
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

public:

	MpSemiImplicitIntegrator() = default;

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

template <int Dim> class MpIntHessianMatrix;

template <int Dim>
class MpIntHessianMatrix : public Eigen::EigenBase<MpIntHessianMatrix<Dim>>
{
	DECLARE_DIM_TYPES(Dim)

public:

	using Scalar = real;
	using RealScalar = real;
	using StorageIndex = int;

	enum
	{
		ColsAtCompileTime = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic,
		IsRowMajor = false
	};

	const GridBasedVectorData<Dim> &velocity;
	const GridBasedVectorData<Dim> &mass;
	const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances;
	const real dt;
	const GridBasedData<Dim, uchar> &collided;

public:

	MpIntHessianMatrix(
		const GridBasedVectorData<Dim> &velocity,
		const GridBasedVectorData<Dim> &mass,
		const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances,
		const real dt,
		const GridBasedData<Dim, uchar> &collided)
		:
		velocity(velocity),
		mass(mass),
		substances(substances),
		dt(dt),
		collided(collided)
	{ }

	virtual ~MpIntHessianMatrix() = default;

	Eigen::Index rows() const { return velocity.count() * Dim; }
	Eigen::Index cols() const { return velocity.count() * Dim; }

	template <typename Rhs>
	Eigen::Product<MpIntHessianMatrix, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs> &x) const
	{
		return Eigen::Product<MpIntHessianMatrix, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
	}
};

template <int Dim>
class MpIntPreconditioner
{
	DECLARE_DIM_TYPES(Dim)

public:

	using StorageIndex = int;

	enum
	{
		ColsAtCompileTime = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic
	};

protected:

	VectorXr _invDiag;

public:

	MpIntPreconditioner() = default;
	virtual ~MpIntPreconditioner() = default;

	explicit MpIntPreconditioner(const MpIntHessianMatrix<Dim> &mat) { compute(mat); }

	Eigen::Index rows() const { return _invDiag.size(); }
	Eigen::Index cols() const { return _invDiag.size(); }

	MpIntPreconditioner &analyzePattern(const MpIntHessianMatrix<Dim> &mat) { return *this; }
	MpIntPreconditioner &compute(const MpIntHessianMatrix<Dim> &mat) { return factorize(mat); }

	MpIntPreconditioner &factorize(const MpIntHessianMatrix<Dim> &mat)
	{
		_invDiag = mat.mass.asVectorXr();

		const real coeff = 4 * mat.dt * mat.velocity.invSpacing() * mat.velocity.invSpacing();
		for (const auto &substance : mat.substances) {
			const real scale = substance->particles.mass() / substance->density() * coeff * coeff;
			substance->particles.forEach([&](const int p) {
				const auto &pos = substance->particles.positions[p];
				const auto dataPoints = mat.velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos);

				std::array<VectorDr, dataPoints.size()> deltaPositions;
				std::array<int, dataPoints.size()> indices;
				// Cache nodes informations.
				for (int i = 0; i < dataPoints.size(); i++) {
					const auto [node, weight] = dataPoints[i];
					deltaPositions[i] = mat.velocity.position(node) - pos;
					indices[i] = int(mat.velocity.index(node)) * Dim;
				}
				// Perform a G2P process.
				MatrixDr weightSum = MatrixDr::Zero();
				for (int i = 0; i < dataPoints.size(); i++) {
					const auto [node, weight] = dataPoints[i];
					if (mat.collided[node]) continue;
					weightSum += weight * VectorDr::Ones() * deltaPositions[i].transpose();
				}
				// Get delta nominal stress tensor.
				const auto deltaStress = substance->computeDeltaStressTensorAtRef(p, weightSum);
				// Perform a P2G process.
				for (int i = 0; i < dataPoints.size(); i++) {
					const auto [node, weight] = dataPoints[i];
					if (mat.collided[node]) continue;
					_invDiag.segment<Dim>(indices[i]) += weight * deltaStress * deltaPositions[i] * scale;
				}
			});
		}

		_invDiag = _invDiag.unaryExpr([](const real x) { return x ? 1 / x : 1; });
		return *this;
	}

	template <typename Rhs, typename Dest>
	void _solve_impl(const Rhs &b, Dest &x) const { x = _invDiag.array() * b.array(); }

	template <typename Rhs>
	inline const Eigen::Solve<MpIntPreconditioner, Rhs> solve(const Eigen::MatrixBase<Rhs> &b) const { return Eigen::Solve<MpIntPreconditioner, Rhs>(*this, b.derived()); }

	Eigen::ComputationInfo info() { return Eigen::Success; }
};

}

namespace Eigen::internal {

template <int Dim> struct traits<PhysX::MpIntHessianMatrix<Dim>> : public Eigen::internal::traits<PhysX::SparseMatrixr> { };

template <int Dim, typename Rhs>
struct generic_product_impl<PhysX::MpIntHessianMatrix<Dim>, Rhs, SparseShape, DenseShape, GemvProduct>
	: generic_product_impl_base<PhysX::MpIntHessianMatrix<Dim>, Rhs, generic_product_impl<PhysX::MpIntHessianMatrix<Dim>, Rhs>>
{
	using Scalar = typename Product<PhysX::MpIntHessianMatrix<Dim>, Rhs>::Scalar;

	template <typename Dest>
	static void scaleAndAddTo(Dest &dst, const PhysX::MpIntHessianMatrix<Dim> &lhs, const Rhs &rhs, const Scalar &alpha)
	{
		dst = lhs.mass.asVectorXr().cwiseProduct(rhs);

		const Scalar coeff = 4 * lhs.dt * lhs.velocity.invSpacing() * lhs.velocity.invSpacing();
		for (const auto &substance : lhs.substances) {
			const Scalar scale = substance->particles.mass() / substance->density() * coeff * coeff;
			substance->particles.forEach([&](const int p) {
				const auto &pos = substance->particles.positions[p];
				const auto dataPoints = lhs.velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos);

				std::array<Matrix<Scalar, Dim, 1>, dataPoints.size()> deltaPositions;
				std::array<int, dataPoints.size()> indices;
				// Cache nodes informations.
				for (int i = 0; i < dataPoints.size(); i++) {
					const auto [node, weight] = dataPoints[i];
					deltaPositions[i] = lhs.velocity.position(node) - pos;
					indices[i] = int(lhs.velocity.index(node)) * Dim;
				}
				// Perform a G2P process.
				Matrix<Scalar, Dim, Dim> weightSum = Matrix<Scalar, Dim, Dim>::Zero();
				for (int i = 0; i < dataPoints.size(); i++) {
					const auto [node, weight] = dataPoints[i];
					if (lhs.collided[node]) continue;
					weightSum += weight * rhs.segment<Dim>(indices[i]) * deltaPositions[i].transpose();
				}
				// Get delta nominal stress tensor.
				const auto deltaStress = substance->computeDeltaStressTensor(p, weightSum);
				// Perform a P2G process.
				for (int i = 0; i < dataPoints.size(); i++) {
					const auto [node, weight] = dataPoints[i];
					if (lhs.collided[node]) continue;
					dst.segment<Dim>(indices[i]) += weight * deltaStress * deltaPositions[i] * scale;
				}
			});
		}
	}
};

}
