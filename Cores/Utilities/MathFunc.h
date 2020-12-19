#pragma once

#include "Utilities/Types.h"

#include <numbers>

#include <cmath>

namespace PhysX::MathFunc
{

template <typename Scalar>
inline constexpr Scalar pow(const Scalar base, const uint exp) { return exp ? base * pow(base, exp - 1) : 1; }

template <typename Scalar>
inline constexpr Scalar square(const Scalar x) { return x * x; }

inline constexpr real dirac(const real x, const real eps)
{
	if (x <= -eps || x >= eps) return 0;
	else return (1 + std::cos(real(std::numbers::pi) * x / eps)) / (2 * eps);
}

inline constexpr real heaviside(const real x, const real eps)
{
	if (x <= -eps) return 0;
	else if (x < eps) return (1 + x / eps + real(std::numbers::inv_pi) * std::sin(real(std::numbers::pi) * x / eps)) / 2;
	else return 1;
}

}
