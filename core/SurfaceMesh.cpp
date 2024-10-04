#include "SurfaceMesh.h"

namespace Pivot {
	void SurfaceMesh::Clear() {
		Positions.clear();
		Normals.clear();
		Indices.clear();
	}

	void SurfaceMesh::Export(std::ostream &out) const {
		IO::Write(out, static_cast<std::uint32_t>(Positions.size()));
		for (auto const &pos : Positions) {
			IO::Write(out, pos.cast<float>().eval());
		}
		IO::Write(out, static_cast<std::uint32_t>(Indices.size()));
		IO::Write(out, Indices);
	}
}
