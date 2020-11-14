#include "StaggeredGrid.h"

namespace PhysX {

template <>
StaggeredGrid<2>::StaggeredGrid(const real spacing, const VectorDi &resolution, const VectorDr &center) :
	_spacing(spacing),
	_resolution(resolution),
	_origin(center - _resolution.cast<real>() * _spacing * real(0.5)),
	_nodeGrid(_spacing, _resolution + VectorDi::Ones(), _origin),
	_cellGrid(_spacing, _resolution, _origin + VectorDr::Ones() * real(0.5) * _spacing),
	_faceGrids {
		Grid<2>(_spacing, _resolution + VectorDi(1, 0), _origin + VectorDr(0, 1) * real(0.5) * _spacing),
		Grid<2>(_spacing, _resolution + VectorDi(0, 1), _origin + VectorDr(1, 0) * real(0.5) * _spacing)
	}
{ }

template <>
StaggeredGrid<3>::StaggeredGrid(const real spacing, const VectorDi &resolution, const VectorDr &center) :
	_spacing(spacing),
	_resolution(resolution),
	_origin(center - _resolution.cast<real>() * _spacing * real(0.5)),
	_nodeGrid(_spacing, _resolution + VectorDi::Ones(), _origin),
	_cellGrid(_spacing, _resolution, _origin + VectorDr::Ones() * real(0.5) * _spacing),
	_faceGrids {
		Grid<3>(_spacing, _resolution + VectorDi(1, 0, 0), _origin + VectorDr(0, 1, 1) * real(0.5) * _spacing),
		Grid<3>(_spacing, _resolution + VectorDi(0, 1, 0), _origin + VectorDr(1, 0, 1) * real(0.5) * _spacing),
		Grid<3>(_spacing, _resolution + VectorDi(0, 0, 1), _origin + VectorDr(1, 1, 0) * real(0.5) * _spacing)
	}
{ }

template class StaggeredGrid<2>;
template class StaggeredGrid<3>;

}
