#include "VolumeSampler.h"

#include "PoissonDiskSampler.h"

namespace Pivot {
	std::vector<Vector3d> VolumeSampler::Sample(Surface const &surface, bool if_Poission) {
		std::vector<Vector3d> samples;

		if(if_Poission){
			PoissonDiskSampler sampler(radius, surface);
			samples = sampler.generatePoints(surface);
		}
		else{
			auto const [minCorner, maxCorner] = surface.GetCornersOfAABB();
			auto const minCell = CellOf(minCorner);
			auto const maxCell = CellOf(maxCorner);
			for (auto const &cell : IterateOverCells(minCell, maxCell)) {
				Vector3d const center = CenterOf(cell);
				Vector3d corner = center;
				for (int axis = 0; axis < 3; ++axis) {
					corner[axis] = corner[axis] - m_CellSize[axis] * .5;
					Vector3d const faceCenter0 = center - Vector3d::Unit(axis) * m_CellSize[axis] * .5;
					if (surface.Surrounds(faceCenter0)) samples.push_back(faceCenter0);
					// if (cell[axis] == maxCell[axis]) {
					// 	Vector2d const faceCenter1 = center + Vector2d::Unit(axis) * m_CellSize[axis] * .5;
					// 	if (surface.Surrounds(faceCenter1)) samples.push_back(faceCenter1);
					// }
				}
				if (surface.Surrounds(corner)) samples.push_back(corner);
			}
		}
		std::cout<<samples.size()<<std::endl;
		return samples;
	}

	std::vector<Vector3i> VolumeSampler::IterateOverCells(Vector3i const &minCell, Vector3i const &maxCell) const {
		std::vector<Vector3i> cells;
		cells.reserve((maxCell.x() - minCell.x() + 2) * (maxCell.y() - minCell.y() + 2) * (maxCell.z() - minCell.z() + 2));
		for (int x = minCell.x(); x <= maxCell.x() + 1; ++x) {
			for (int y = minCell.y(); y <= maxCell.y() + 1; ++y) {
				for (int z = minCell.z(); z <= maxCell.z() + 1; ++z)
					cells.emplace_back(x, y, z);
			}
		}
		return cells;
	}
}