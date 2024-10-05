#pragma once

#include "Surface.h"

namespace Pivot {
	class SurfaceMesh : public Surface {
	public:
		SurfaceMesh() = default;

		virtual Vector3d ClosestNormalOf (Vector3d const &pos) const override { return Vector3d::Zero(); } // FIXME
		virtual double   SignedDistanceTo(Vector3d const &pos) const override { return 0; } // FIXME
        virtual std::pair<Vector3d, Vector3d> GetCornersOfAABB() const override { return {Vector3d::Zero(), Vector3d::Zero()}; }

		void Clear();
		void Export(std::ostream &out) const;

		void ComputeNormals();

	public:
		std::vector<Vector3d>      Positions;
		std::vector<Vector3d>      Normals;
		std::vector<std::uint32_t> Indices;
	};
}
