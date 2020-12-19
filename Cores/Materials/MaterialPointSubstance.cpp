#include "MaterialPointSubstance.h"

namespace PhysX {

template <int Dim>
void MaterialPointSubstance<Dim>::write(std::ofstream &fout) const
{
	IO::writeValue(fout, uint(particles.size()));
	particles.forEach([&](const int i) {
		IO::writeValue(fout, particles.positions[i].template cast<float>().eval());
	});
	if constexpr (Dim == 3) {
		particles.forEach([&](const int i) {
			IO::writeValue(fout, VectorDf::Unit(2).eval());
		});
	}
}

template <int Dim>
void MaterialPointSubstance<Dim>::save(std::ofstream &fout) const
{
	IO::writeValue(fout, uint(particles.size()));
	particles.positions.save(fout);
}

template <int Dim>
void MaterialPointSubstance<Dim>::load(std::ifstream &fin)
{
	uint particlesCnt;
	IO::readValue(fin, particlesCnt);
	particles.resize(particlesCnt);
	particles.positions.load(fin);
	reinitialize();
}

template <int Dim>
void MaterialPointSubstance<Dim>::reinitialize()
{
	velocities.resize(&particles);
	velocities.setZero();

	velocityDerivatives.resize(&particles);
	velocityDerivatives.setZero();
}


template class MaterialPointSubstance<2>;
template class MaterialPointSubstance<3>;

}
