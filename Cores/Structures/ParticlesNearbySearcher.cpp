#include "ParticlesNearbySearcher.h"

#include "Utilities/MathFunc.h"

#include <algorithm>
#include <array>

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
	int keysCnt = 0;
	std::array<int, MathFunc::pow(3, Dim)> hashKeys;
	for (const auto &coord : _grid.quadraticNearbyDataPoints(pos))
		hashKeys[keysCnt++] = getHashKey(coord);
	std::sort(hashKeys.begin(), hashKeys.end());

	for (int i = 0; i < keysCnt; i++) {
		if (i && hashKeys[i] == hashKeys[size_t(i) - 1]) continue;
		for (const int j : _buckets[hashKeys[i]]) {
			if ((positions[j] - pos).squaredNorm() < _squaredRadius) func(j, positions[j]);
		}
	}
}

template class ParticlesNearbySearcher<2>;
template class ParticlesNearbySearcher<3>;

template class HashGridSearcher<2>;
template class HashGridSearcher<3>;

}
