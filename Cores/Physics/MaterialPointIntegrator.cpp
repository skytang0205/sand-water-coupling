#include "MaterialPointIntegrator.h"

namespace PhysX {

template <int Dim>
void MpSymplecticEulerIntegrator<Dim>::integrate(
	GridBasedVectorData<Dim> &momentum,
	const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances,
	const real dt)
{
	for (auto &substance : _substances) {
		const real stressCoeff = -dt * 4 * momentum.invSpacing() * momentum.invSpacing() * substance->particles.mass() / substance->density();
		substance->computeStressTensors(stresses);

		substance->particles.forEach([&](const int i) {
			const VectorDr pos = substance->particles.positions[i];
			const MatrixDr stress = stresses[i] * stressCoeff;
			for (const auto [node, weight] : momentum.grid()->quadraticBasisSplineIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity.position(node) - pos;
				momentum[node] += stress * deltaPos * weight;
			}
		});
	}
}

}
