#include "VectorGridField.h"

namespace PhysX {

template <int Dim>
Vector<real, Dim> VectorGridField<Dim, FaceCentered>::operator()(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;

	VectorDr vec = VectorDr::Zero();
	for (int axis = 0; axis < Dim; axis++) {
		getFaceLerp(axis, pos, coords, weights);
		for (int i = 0; i < (1 << Dim); i++)
			vec[axis] += _data[axis][coords[i]] * weights[i];
	}
	return vec;
}

template <int Dim>
real VectorGridField<Dim, FaceCentered>::divergenceAtCellCenter(const VectorDi &cell) const
{
	real val = 0;
	for (int i = 0; i < Dim; i++) {
		val += _data[i][cell + VectorDi::Unit(i)] - _data[i][cell];
	}
	return val / (2 * _spacing);
}

template <int Dim>
real VectorGridField<Dim, FaceCentered>::divergence(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	getCellLerp(pos, coords, weights);

	real val = 0;
	for (int i = 0; i < (1 << Dim); i++)
		val += divergenceAtCellCenter(coords[i]) * weights[i];
	return val;
}

template class VectorGridField<2, FaceCentered>;
template class VectorGridField<3, FaceCentered>;

}
