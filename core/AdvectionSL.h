#pragma once

#include "TriLerp.h"
#include "SGridData.h"

namespace Pivot {
	class AdvectionSL {
	public:
		template <int RkOrder, typename Type>
			requires (1 <= RkOrder && RkOrder <= 4)
		static void Solve(GridData<Type> &grData, SGridData<double> const &flow, double dt) {
			GridData<Type> newGrData(grData.GetGrid());
			ParallelForEach(grData.GetGrid(), [&](Vector3i const &coord) {
				Vector3d const pos = grData.GetGrid().PositionOf(coord);
				newGrData[coord] = TriLerp::Interpolate(grData, Trace<RkOrder>(pos, flow, -dt));
			});
			grData = newGrData;
		}

		template <int RkOrder, typename Type>
			requires (1 <= RkOrder && RkOrder <= 4)
		static void Solve(SGridData<Type> &sgrData, SGridData<double> const &flow, double dt) {
			SGridData<Type> newSgrData(sgrData.GetGrids());
			ParallelForEach(sgrData.GetGrids(), [&](int axis, Vector3i const &face) {
				Vector3d const pos = sgrData[axis].GetGrid().PositionOf(face);
				newSgrData[axis][face] = TriLerp::Interpolate(sgrData[axis], Trace<RkOrder>(pos, flow, -dt));
			});
			sgrData = newSgrData;
		}

		template <int RkOrder>
			requires (1 <= RkOrder && RkOrder <= 4)
		static Vector3d Trace(Vector3d const &startPos, SGridData<double> const &flow, double dt) {
			if constexpr (RkOrder == 1) {
				return startPos + TriLerp::Interpolate(flow, startPos) * dt;
			} else if constexpr (RkOrder == 2) { // the midpoint scheme
				Vector3d const vel0 = TriLerp::Interpolate(flow, startPos);
				Vector3d const vel1 = TriLerp::Interpolate(flow, startPos + vel0 * dt * .5);
				return startPos + vel1 * dt;
			} else if constexpr (RkOrder == 3) { // the classical third order Runge-Kutta scheme
				Vector3d const vel0 = TriLerp::Interpolate(flow, startPos);
				Vector3d const vel1 = TriLerp::Interpolate(flow, startPos + vel0 * dt * .5);
				Vector3d const vel2 = TriLerp::Interpolate(flow, startPos + (2 * vel1 - vel0) * dt);
				return startPos + (vel0 + 4 * vel1 + vel2) * dt / 6;
			} else if constexpr (RkOrder == 4) { // the classical fourth order Runge-Kutta scheme
				Vector3d const vel0 = TriLerp::Interpolate(flow, startPos);
				Vector3d const vel1 = TriLerp::Interpolate(flow, startPos + vel0 * dt * .5);
				Vector3d const vel2 = TriLerp::Interpolate(flow, startPos + vel1 * dt * .5);
				Vector3d const vel3 = TriLerp::Interpolate(flow, startPos + vel2 * dt);
				return startPos + (vel0 + 2 * vel1 + 2 * vel2 + vel3) * dt / 6;
			}
		}
	};
}
