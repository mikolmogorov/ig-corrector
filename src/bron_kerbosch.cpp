#include "bron_kerbosch.h"
#include <iostream>
#include <algorithm>

namespace
{
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

	void bron_kerbosch_pivot(const Graph::gr_size_type& threshold, const Graph& g, 
							std::vector<Graph::vertex_set>& vs,
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
			auto u = g.pivot(p, x);

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
}

Graph::vert_id Graph::pivot(const Graph::vertex_set& p, const Graph::vertex_set& x) const
{
	Graph::vertex_set base;
	for (auto it = p.begin(); it != p.end(); ++it)
		base.insert(*it);
	for (auto it = x.begin(); it != x.end(); ++it)
		base.insert(*it);

	Graph::vert_id picked = Graph::null_vertex();
	for (auto it = base.begin(); it != base.end(); ++it) 
	{
		if (picked == Graph::null_vertex() || this->degree_size(picked) < this->degree_size(*it)) 
		{
	 		picked = *it;
		}
	}
	return picked;
}

void Graph::find_cliques(std::vector<Graph::vertex_set>& vs) const
{
	vertex_set r, p, x;
	for(auto &v : *this)
	{
		p.insert(v.first);
	}
	bron_kerbosch_pivot(2, *this, vs, r, p, x);
}
