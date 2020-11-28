#pragma once

#include "Utilities/Types.h"

#include <numbers>

#include <cmath>

namespace PhysX::MathFunc
{

inline constexpr real cubicCatmullRomSplineCoefficient(const int i, const real x)
{
	switch (i) {
	case 0: return -x / 3 + x * x / 2 - x * x * x / 6;
	case 1: return 1 - x * x + (x * x * x - x) / 2;
	case 2: return x + (x * x - x * x * x) / 2;
	case 3: return (x * x * x - x) / 6;
	default: return 0;
	}
}

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
