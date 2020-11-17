#pragma once

#include "Geometries/GridBasedImplicitSurface.h"

namespace PhysX {

template <int Dim>
class LevelSetReinitializer
{
	DECLARE_DIM_TYPES(Dim)

protected:

	int _reinitMaxIters;

public:

	LevelSetReinitializer(const int reinitMaxIters) : _reinitMaxIters(reinitMaxIters) { }

	LevelSetReinitializer(const LevelSetReinitializer &rhs) = delete;
	LevelSetReinitializer &operator=(const LevelSetReinitializer &rhs) = delete;
	virtual ~LevelSetReinitializer() = default;

	virtual void reinitialize(GridBasedImplicitSurface<Dim> &levelSet) = 0;
};

template <int Dim>
class FastMarchingReinitializer final : public LevelSetReinitializer<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using LevelSetReinitializer<Dim>::_reinitMaxIters;

	GridBasedData<Dim> _tent; // tentative signed distance
	GridBasedData<Dim, uchar> _valid;

public:

	FastMarchingReinitializer(const Grid<Dim> *const grid, const int reinitMaxIters);

	virtual void reinitialize(GridBasedImplicitSurface<Dim> &levelSet) override;

protected:

	void initInterface(const GridBasedScalarField<Dim> &phi);
	void fastMarching();
};

}
