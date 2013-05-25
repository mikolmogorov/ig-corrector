#ifndef BRON_KERBOSH_H
#define BRON_KERBOSH_H

#include <vector>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

/*
	A pivoting version of Bron-Kerbosch algorithm
	
	Adapted from: git.io/bk73

## References:

	- Bron, C. and Kerbosch, J. 1973.
	Algorithm 457: finding all cliques of an undirected graph.
	(http://dl.acm.org/citation.cfm?id=362367)
	- Etsuji Tomita, Akira Tanaka, Haruhisa Takahashi,
	The worst-case time complexity for generating all maximal
	cliques and computational experiments.
	(http://dl.acm.org/citation.cfm?id=1217600)
	- F. Cazals, C. Karande,
	A note on the problem of reporting maximal cliques.
	(http://dl.acm.org/citation.cfm?id=1450505)

*/

class Graph 
{
public:
	typedef size_t gr_size_type;
	typedef size_t vert_id;
	typedef std::ptrdiff_t vertex_difference;
	typedef std::unordered_set<vert_id> vertex_set;
	typedef std::unordered_map<vert_id, vertex_set>::const_iterator vertex_set_iter;

	Graph (): nbrs_()
	{}

	vertex_set_iter begin() const
	{
		return nbrs_.begin();
	}

	vertex_set_iter end() const
	{
		return nbrs_.end();
	}

	gr_size_type size() const
	{
		return nbrs_.size();
	}

	static vert_id null_vertex()
	{
		return std::numeric_limits<vert_id>::max();
	}

	gr_size_type degree_size(const vert_id& v) const
	{
		auto it = nbrs_.find(v);
		if (it == nbrs_.end()) return false;
		return it->second.size();
	}

	bool is_connected(const vert_id& v1, const vert_id& v2) const
	{
		auto it = nbrs_.find(v1);
		if (it == nbrs_.end()) return false;
		return (it->second.find(v2) != it->second.end());
	}
/*
	void load(std::istream& in)
	{
		std::string line;
		while (getline(in, line)) 
		{
			if (line.empty()) continue;

			std::stringstream ss(line);
			vert_id cost, src, dst;
			ss >> cost >> src >> dst;
			add_edge(src, dst);
		}
	}
*/
	void add_edge(const vert_id& src, const vert_id& dst)
	{
		nbrs_[src].insert(dst);
	}

	void find_cliques(std::vector<Graph::vertex_set>& vs) const;

	vert_id pivot(const Graph::vertex_set& p, const Graph::vertex_set& x) const;
private:
	std::unordered_map<vert_id, std::unordered_set<vert_id>> nbrs_;
}; // class Graph


#endif
