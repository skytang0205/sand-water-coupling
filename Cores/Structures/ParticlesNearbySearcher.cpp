#include "ParticlesNearbySearcher.h"

namespace PhysX {

template <int Dim>
void HashGridSearcher<Dim>::reset(const ParticlesVectorAttribute<Dim> &positions)
{
	_buckets.resize(positions.size() * 2);
	for (auto &bucket : _buckets) bucket.clear();
	positions.forEach([&](const int i) {
		_buckets[getHashKey(positions[i])].push_back(i);
	});
}

template <int Dim>
void HashGridSearcher<Dim>::forEach(const ParticlesVectorAttribute<Dim> &positions, const VectorDr &pos, const std::function<void(const int, const VectorDr &)> &func)
{
	for (const auto &coord : _grid.quadraticNearbyDataPoints(pos)) {
		for (const int j : _buckets[getHashKey(coord)]) {
			if ((positions[j] - pos).squaredNorm() < _squaredRadius) func(j, positions[j]);
		}
	}
}

template class ParticlesNearbySearcher<2>;
template class ParticlesNearbySearcher<3>;

template class HashGridSearcher<2>;
template class HashGridSearcher<3>;

}
