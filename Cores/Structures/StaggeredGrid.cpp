#include "StaggeredGrid.h"

namespace PhysX {

template <>
StaggeredGrid<2>::StaggeredGrid(const real spacing, const VectorDi &resolution, const VectorDr &center) :
	_spacing(spacing),
	_invSpacing(1 / _spacing),
	_resolution(resolution + VectorDi::Ones() * _kBoundaryWidth * 2),
	_origin(center - _resolution.template cast<real>() * _spacing / 2),
	_nodeGrid(_spacing, _resolution + VectorDi::Ones(), _origin),
	_cellGrid(_spacing, _resolution, _origin + VectorDr::Ones() * _spacing / 2),
	_faceGrids {
		Grid<2>(_spacing, _resolution + VectorDi(1, 0), _origin + VectorDr(0, 1) * _spacing / 2),
		Grid<2>(_spacing, _resolution + VectorDi(0, 1), _origin + VectorDr(1, 0) * _spacing / 2)
	}
{ }

template <>
StaggeredGrid<3>::StaggeredGrid(const real spacing, const VectorDi &resolution, const VectorDr &center) :
	_spacing(spacing),
	_invSpacing(1 / _spacing),
	_resolution(resolution + VectorDi::Ones() * _kBoundaryWidth * 2),
	_origin(center - _resolution.template cast<real>() * _spacing / 2),
	_nodeGrid(_spacing, _resolution + VectorDi::Ones(), _origin),
	_cellGrid(_spacing, _resolution, _origin + VectorDr::Ones() * _spacing / 2),
	_faceGrids {
		Grid<3>(_spacing, _resolution + VectorDi(1, 0, 0), _origin + VectorDr(0, 1, 1) * _spacing / 2),
		Grid<3>(_spacing, _resolution + VectorDi(0, 1, 0), _origin + VectorDr(1, 0, 1) * _spacing / 2),
		Grid<3>(_spacing, _resolution + VectorDi(0, 0, 1), _origin + VectorDr(1, 1, 0) * _spacing / 2)
	}
{ }

template class StaggeredGrid<2>;
template class StaggeredGrid<3>;

}
