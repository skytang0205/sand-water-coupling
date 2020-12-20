#pragma once

#include "Materials/FixedCorotatedModel.h"
#include "Materials/NeoHookeanModel.h"

#include <concepts>

namespace PhysX {

template <typename Model, int Dim>
concept Hyperelastic = requires (const Matrix<Dim, real> &defmGrad, const real lambda, const real mu) {
	{ Model::computeEnergy(defmGrad, lambda, mu) } -> std::convertible_to<real>;
	{ Model::computeNominalStressTensor(defmGrad, lambda, mu) } -> std::convertible_to<Matrix<Dim, real>>;
	{ Model::computeStressTensor(defmGrad, lambda, mu) } -> std::convertible_to<Matrix<Dim, real>>;
	{ Model::computeStressTensorMultipliedByJ(defmGrad, lambda, mu) } -> std::convertible_to<Matrix<Dim, real>>;
	{ Model::computeEnergyHessian(defmGrad, lambda, mu) } -> std::convertible_to<Matrix<Dim * Dim, real>>;
};

}
