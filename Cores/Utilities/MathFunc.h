#pragma once

#include "Utilities/Types.h"

#include <numbers>

#include <cmath>

namespace PhysX::MathFunc
{

inline constexpr int pow(int base, uint exp) { return exp ? base * pow(base, exp - 1) : 1; }

inline constexpr real quadraticBasisSplineCoefficient(real x)
{
	if (x < 0) x = -x;
	if (x < real(1) / 2) return real(3) / 4 - x * x;
	else if (x < real(3) / 2) return (real(3) / 2 - x) * (real(3) / 2 - x) / 2;
	else return 0;
}

inline constexpr real cubicBasisSplineCoefficient(real x)
{
	if (x < 0) x = -x;
	if (x < 1) return x * x * x / 2 - x * x + real(2) / 3;
	else if (x < 2) return (2 - x) * (2 - x) * (2 - x) / 6;
	else return 0;
}

inline constexpr real cubicCatmullRomCoefficient(const int i, const real x)
{
	// TODO: fix range of x
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
