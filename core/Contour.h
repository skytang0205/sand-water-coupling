#pragma once

#include "GridData.h"
#include "SurfaceMesh.h"

namespace Pivot {
	class Contour {
	public:
		Contour(Grid const &nodeGrid);

		SurfaceMesh const &GetMesh() const { return m_Mesh; }

		void Generate(GridData<double> const &grData, double value = 0.);

		int  VertexIndexOf(int axis, Vector2i const &edge) { return m_EdgeMark[axis][edge]; }

	private:
		#include "MarchingCubesTables.inc"

		Grid                         const &m_NodeGrid;
		Grid                                m_CellGrid;
		std::array<Grid, 2>                 m_EdgeGrids;
		std::array<GridData<int>, 2>        m_EdgeMark;

		SurfaceMesh                         m_Mesh;
	};
}
