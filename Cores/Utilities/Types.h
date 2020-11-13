#pragma once

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

#include <cstddef>

namespace PhysX {

#define DECLARE_EIGEN_VECTOR_TYPES(type, t)							\
using	Vector2##t			=	Eigen::Vector2##t;					\
using	Vector3##t			=	Eigen::Vector3##t;					\
using	Vector4##t			=	Eigen::Vector4##t;					\
using	VectorX##t			=	Eigen::VectorX##t;

#define DECLARE_EIGEN_MATRIX_TYPES(type, t)							\
using	Matrix2##t			=	Eigen::Matrix2##t;					\
using	Matrix3##t			=	Eigen::Matrix3##t;					\
using	Matrix4##t			=	Eigen::Matrix4##t;					\
using	SparseMatrix##t		=	Eigen::SparseMatrix<type>;

DECLARE_EIGEN_VECTOR_TYPES(int, i)
DECLARE_EIGEN_VECTOR_TYPES(float, f)
DECLARE_EIGEN_VECTOR_TYPES(double, d)

DECLARE_EIGEN_MATRIX_TYPES(float, f)
DECLARE_EIGEN_MATRIX_TYPES(double, d)

#undef DECLARE_EIGEN_VECTOR_TYPES
#undef DECLARE_EIGEN_MATRIX_TYPES

#define DECLARE_REAL_TYPES(type, t)									\
using	real				=	type;								\
using	Vector2r			=	Vector2##t;							\
using	Vector3r			=	Vector3##t;							\
using	Vector4r			=	Vector4##t;							\
using	VectorXr			=	VectorX##t;							\
using	Matrix2r			=	Matrix2##t;							\
using	Matrix3r			=	Matrix3##t;							\
using	SparseMatrixr		=	SparseMatrix##t;

#ifdef USE_FLOAT
DECLARE_REAL_TYPES(float, f)
#else
DECLARE_REAL_TYPES(double, d)
#endif

#undef DECLARE_REAL_TYPES

template <class Scalar, int Dim> using Vector = Eigen::Matrix<Scalar, Dim, 1>;
template <class Scalar, int Dim> using Matrix = Eigen::Matrix<Scalar, Dim, Dim>;

#define DECLARE_DIM_TYPES(Dim)										\
static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");	\
using	VectorDr			=	Vector<real, Dim>;					\
using	VectorDf			=	Vector<float, Dim>;					\
using	VectorDi			=	Vector<int, Dim>;					\
using	MatrixDr			=	Matrix<real, Dim>;

using	uchar				=	unsigned char;
using	ushort				=	unsigned short;
using	uint				=	unsigned int;
using	llong				=	long long;
using	ullong				=	unsigned long long;

using	std::size_t;

} // namespaxe PhysX
