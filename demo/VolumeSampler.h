#pragma once

#include "Surface.h"

namespace Pivot {
	class VolumeSampler {
	public:
		VolumeSampler(double latticeConstant) : radius(latticeConstant), m_CellSize { latticeConstant * std::numbers::sqrt2, latticeConstant * std::numbers::sqrt2, latticeConstant * std::numbers::sqrt2 } { }

		std::vector<Vector3d> Sample(Surface const &surface, bool if_Poission);
	
	private:
		Vector3i CellOf  (Vector3d const &pos ) const { return (pos.array() / m_CellSize - c_Offset).floor().template cast<int>().matrix(); }
		Vector3d CenterOf(Vector3i const &cell) const { return (cell.template cast<double>() + Vector3d::Constant(.5 + c_Offset)).cwiseProduct(m_CellSize.matrix()); }

		std::vector<Vector3i> IterateOverCells(Vector3i const &minIndex, Vector3i const &maxIndex) const;
	
	private:
		static constexpr double c_Offset = .02407112028; // A magic number to avoid conflicts

		double const radius;

		Array3d const m_CellSize;
	};
}