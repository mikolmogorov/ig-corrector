
#include "hierarchial_impl.h"
#include "edit_distance.h"
#include <algorithm>
#include <iostream>
#include <utility>

namespace
{
	template <class T>
	std::pair<T, T> ordered_pair(const T& e1, const T& e2)
	{
		if (e1 <= e2)
		{
			return std::make_pair(e1, e2);
		}
		else
		{
			return std::make_pair(e2, e1);
		}
	}
}

void HierarchialClust::cluster(const FastaSequences& seqs, FastaSet& output, float cutoff)
{
	_cutoff = cutoff;

	int counter = 0;
	for (auto seq : seqs)
	{
		_fastaHash[seq.description] = seq.sequence;
		_clusters[counter].push_back(seq.description);
		counter += 1;
	}

	for (int i = 0; i < counter; ++i)
	{
		for (int j = i + 1; j < counter; ++j)
		{
			float dist = clusterDist(i, j);
			_distances[ordered_pair(i, j)] = dist;
		}
	}

	while (_clusters.size() > 1)
	{
		bool ret = this->step();
		if (!ret) break;
	}

	this->outputClusters(output);
}

bool HierarchialClust::step()
{
	auto compare = [](const DistHash::value_type& p1, const DistHash::value_type& p2)
		{return p1.second < p2.second;};
	KeyPair minKey = std::min_element(_distances.begin(), _distances.end(), compare)->first;
	if (_distances[minKey] > _cutoff) return false;

	std::cerr << "merging " << minKey.first << " " << minKey.second << " " << _distances[minKey] << std::endl;
	//merge two clusters
	auto itList = _clusters.find(minKey.first);
	itList->second.splice(itList->second.end(), _clusters[minKey.second]);

	//update distances
	for (auto clust : _clusters)
	{
		if (clust.first != minKey.second)
		{
			_distances.erase(ordered_pair(clust.first, minKey.second));
		}
		if (clust.first != minKey.first && clust.first != minKey.second)
		{
			_distances[ordered_pair(clust.first, minKey.first)] = this->clusterDist(clust.first, minKey.first);
		}
	}
	_clusters.erase(minKey.second);
}

int HierarchialClust::sequenceDistance(const std::string& seq1, const std::string& seq2)
{
	static std::unordered_map<std::pair<std::string, std::string>, int, pair_hash<std::string>> distCache;

	auto key = ordered_pair(seq1, seq2);
	auto itFound = distCache.find(key);
	if (itFound == distCache.end())
	{
		int distance = edit_distance(_fastaHash[seq1], _fastaHash[seq2]);
		distCache[key] = distance;
		return distance;
	}
	return itFound->second;
}

float HierarchialClust::clusterDist(int clust1, int clust2)
{
	float dist = 0.0f;
	for (std::string& seq1 : _clusters[clust1])
	{
		for (std::string& seq2 : _clusters[clust2])
		{
			dist += this->sequenceDistance(seq1, seq2);
		}
	}
	return dist / (_clusters[clust1].size() * _clusters[clust2].size());
}

void HierarchialClust::outputClusters(FastaSet& out)
{
	int id = 0;
	for (auto clust : _clusters)
	{
		out.push_back(std::vector<FastaRecord>());
		for (auto head : clust.second)
		{
			out.back().push_back(FastaRecord(_fastaHash[head], head, id++));
		}
	}
}
