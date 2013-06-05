#ifndef COMMON_H
#define COMMON_H

#include <utility>
#include <unordered_map>

template <class T>
struct pair_hash 
{
	size_t operator() (const std::pair<T, T>& v) const 
	{
		return std::hash<T>()(v.first) * 31 + std::hash<T>()(v.second);
	}
};

#endif
