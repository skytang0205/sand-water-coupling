#pragma once

#include "GridData.h"

namespace Pivot {
	template <int Order>
	class TriBSpline {
	public:
		static std::array<Vector3i, (Order + 1) * (Order + 1) * (Order + 1)> GetPoints(Grid const &grid, Vector3d const &pos) {
			std::array<Vector3i, (Order + 1) * (Order + 1) * (Order + 1)> pts;
			Vector3i const lower = grid.CalcLower<Order>(pos);
			for (int i = 0; i <= Order; i++) {
				for (int j = 0; j <= Order; j++) {
					for (int k = 0; k <= Order; k++) {
						pts[(i * (Order + 1) + j) * (Order + 1) + k] = lower + Vector3i(i, j, k);
					}
				}
			}
			return pts;
		}
	};
}