
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#define DEBUG

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

	void add_edge(const vert_id& src, const vert_id& dst)
	{
		nbrs_[src].insert(dst);
	}

private:
	std::unordered_map<vert_id, std::unordered_set<vert_id>> nbrs_;
}; // class Graph

Graph::vert_id pivot(const Graph& g, const Graph::vertex_set& p, const Graph::vertex_set& x)
{
	Graph::vertex_set base;
	for (auto it = p.begin(); it != p.end(); ++it)
		base.insert(*it);
	for (auto it = x.begin(); it != x.end(); ++it)
		base.insert(*it);

	Graph::vert_id picked = Graph::null_vertex();
	for (auto it = base.begin(); it != base.end(); ++it) 
	{
		if (picked == Graph::null_vertex() || g.degree_size(picked) < g.degree_size(*it)) 
		{
	 		picked = *it;
		}
	}
	return picked;
}

void nbr_intersection(const Graph& g, Graph::vert_id v, const Graph::vertex_set& s, 
						Graph::vertex_set& ret)
{
	for (auto it = s.begin(); it != s.end(); ++it) 
	{
    	if (g.is_connected(v, *it)) 
		{
      		ret.insert(*it);
    	}
  	}
}

void bron_kerbosch_pivot(const Graph::gr_size_type& threshold, const Graph& g, std::vector<Graph::vertex_set>& vs,
    					Graph::vertex_set& r, Graph::vertex_set& p, Graph::vertex_set& x)
{
	if (p.empty() && x.empty()) 
  	{
    	if (r.size() > threshold) 
		{
      		vs.push_back(Graph::vertex_set());
      		std::copy(r.begin(), r.end(), inserter(vs.back(), vs.back().begin()));
    	}
  	} 
	else 
	{
    	auto u = pivot(g, p, x);

    	Graph::vertex_set rems;
   		for (auto it = p.begin(); it != p.end(); ++it)
      		if (!g.is_connected(u, *it))
        		rems.insert(*it);

    	for (auto it = rems.begin(); it != rems.end(); ++it) 
		{
      		Graph::vert_id v = *it;
      		Graph::vertex_set new_r, new_p, new_x;
      		std::copy(r.begin(), r.end(), inserter(new_r, new_r.begin()));
      		new_r.insert(v);
      		nbr_intersection(g, v, p, new_p);
      		nbr_intersection(g, v, x, new_x);

      		bron_kerbosch_pivot(threshold, g, vs, new_r, new_p, new_x);

      		p.erase(p.find(v));
      		x.insert(v);
    	}
  	}
}

void find_cliques(const Graph& g, std::vector<Graph::vertex_set>& vs)
{
	Graph::vertex_set r, p, x;
	for(auto &v : g)
	{
		p.insert(v.first);
	}
	bron_kerbosch_pivot(2, g, vs, r, p, x);
}

int main()
{
	Graph g;

#ifdef DEBUG
  /* sample_input:
   *   12 - 11   1425
   *   |  x  | /  |
   *   8 -- 10 - 1426
   *
   */
	std::stringstream sample_input;
	sample_input << "1       8       10"     << std::endl
               	<< "1       8       11"     << std::endl
               << "1       8       12"     << std::endl
               << "1       10      8"      << std::endl
               << "1       10      11"     << std::endl
               << "1       10      12"     << std::endl
               << "1       10      1425"   << std::endl
               << "1       10      1426"   << std::endl
               << "1       11      8"      << std::endl
               << "1       11      10"     << std::endl
               << "1       11      12"     << std::endl
               << "1       12      8"      << std::endl
               << "1       12      10"     << std::endl
               << "1       12      11"     << std::endl
               << "1       1425    10"     << std::endl
               << "1       1425    1426"   << std::endl
               << "1       1426    10"     << std::endl
               << "1       1426    1425"   << std::endl;
	g.load(sample_input);
#else
	g.load(std::cin);
#endif

	std::vector<Graph::vertex_set> vvs;
	find_cliques(g, vvs);

	std::for_each(
      vvs.begin(),
      vvs.end(),
      [&](const Graph::vertex_set& vs) {
        std::cout << "clique:";
        std::for_each(
            vs.begin(),
            vs.end(),
            [&](const Graph::vert_id v) {
              std::cout << ' ' << v;
            });
        std::cout << std::endl;
      });
}
