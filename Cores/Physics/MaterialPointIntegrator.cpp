#include "MaterialPointIntegrator.h"

#include "Solvers/IterativeSolver.h"

namespace PhysX {

template <int Dim>
void MpSemiImplicitIntegrator<Dim>::integrate(
	GridBasedVectorData<Dim> &velocity,
	const GridBasedScalarData<Dim> &mass,
	const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances,
	const real dt,
	const GridBasedData<Dim, uchar> &collided)
{
	GridBasedVectorData<Dim> massExt(mass.grid());
	mass.parallelForEach([&](const VectorDi &node) {
		massExt[node] = VectorDr::Constant(mass[node]);
	});

	IterativeSolver::solve<MpIntHessianMatrix<Dim>, Eigen::BiCGSTAB<MpIntHessianMatrix<Dim>, MpIntPreconditioner<Dim>>>(
		MpIntHessianMatrix<Dim>(velocity, massExt, substances, dt, collided),
		velocity.asVectorXr(),
		massExt.asVectorXr().cwiseProduct(velocity.asVectorXr()));
}

template class MaterialPointIntegrator<2>;
template class MaterialPointIntegrator<3>;
template class MpSymplecticEulerIntegrator<2>;
template class MpSymplecticEulerIntegrator<3>;
template class MpSemiImplicitIntegrator<2>;
template class MpSemiImplicitIntegrator<3>;

}
