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
	MaterialPointIntegrator(const MaterialPointIntegrator & rhs) = delete;
	MaterialPointIntegrator &operator=(const SpringMassSysIntegrator & rhs) = delete;
	virtual ~MaterialPointIntegrator() = default;

	virtual void integrate(
		GridBasedVectorData<Dim> &momentum,
		const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances,
		const real dt) = 0;
};

template <int Dim>
class MpSymplecticEulerIntegrator : public MaterialPointIntegrator<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	ParticlesBasedData<Dim, MatrixDr> _stresses;

public:

	MpSymplecticEulerIntegrator() = default;
	MpSymplecticEulerIntegrator(const MpSymplecticEulerIntegrator & rhs) = delete;
	MpSymplecticEulerIntegrator &operator=(const MpSymplecticEulerIntegrator & rhs) = delete;
	virtual ~MpSymplecticEulerIntegrator() = default;

	virtual void integrate(
		GridBasedVectorData<Dim> &momentum,
		const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances,
		const real dt) override;
};

}
