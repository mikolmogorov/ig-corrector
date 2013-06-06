#include <iostream>
#include <limits>
#include <functional>
#include <algorithm>
#include "xalign_impl.h"

namespace
{
	const std::unordered_map<char, std::list<std::string>> backCodonTable = 
	{
		{'A', {"GCT", "GCC", "GCA", "GCG"}},
		{'R', {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}},
		{'N', {"AAT", "AAC"}},
		{'D', {"GAT", "GAC"}},
		{'C', {"TGT", "TGC"}},
		{'Q', {"CAA", "CAG"}},
		{'E', {"GAA", "GAG"}},
		{'G', {"GGT", "GGC", "GGA", "GGG"}},
		{'H', {"CAT", "CAC"}},
		{'I', {"ATT", "ATC", "ATA"}},
		{'L', {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}},
		{'K', {"AAA", "AAG"}},
		{'M', {"ATG"}},
		{'F', {"TTT", "TTC"}},
		{'P', {"CCT", "CCC", "CCA", "CCG"}},
		{'S', {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}},
		{'T', {"ACT", "ACC", "ACA", "ACG"}},
		{'W', {"TGG"}},
		{'Y', {"TAT", "TAC"}},
		{'V', {"GTT", "GTC", "GTA", "GTG"}},
		{'*', {"TAA", "TGA", "TAG"}}
	};
}

void Xaligner::setQuerry(const std::string& querry)
{
	this->clearGraph();
	GraphNode* graph = this->buildGraph(querry);
	_root = graph;
}

void Xaligner::clearGraph()
{
	if (_root == nullptr) return;
	
	std::vector<GraphNode*> graphNodes;
	topSort(_root, graphNodes);
	for (GraphNode* node : graphNodes)
	{
		delete node;
	}
}

void Xaligner::align(const std::string& sequence, int threshold, std::vector<AlignInfo>& result)
{
	if (_root == nullptr) return;
	this->alignOnGraph(sequence, _root, threshold, result);
}

Xaligner::GraphNode* Xaligner::buildGraph(const std::string& sequence)
{
	GraphNode* root = new GraphNode();
	GraphNode* currentRoot = root;

	std::string variantStr;
	bool isVariant = false;
	for (char sym : sequence)
	{
		GraphNode* newRoot = new GraphNode();
		std::list<std::string> codons;
		if (sym == '[')
		{
			isVariant = true;
			variantStr = "";
			continue;
		}
		else if (sym == ']')
		{
			for (char n : variantStr)
			{
				codons.insert(codons.end(), backCodonTable.at(n).begin(), backCodonTable.at(n).end());
			}
			isVariant = false;
		}
		else if (isVariant)
		{
			variantStr += sym;
			continue;
		}
		else
		{
			codons = backCodonTable.at(sym);
		}
		
		for (std::string& codon : codons)
		{
			GraphNode* node = currentRoot;
			for (size_t i = 0; i < codon.length(); ++i)
			{
				char sym = codon[i];
				if (i < 2)
				{
					auto itEdges = node->outEdges.find(sym);
					if (itEdges == node->outEdges.end())
					{
						node->outEdges[sym] = new GraphNode();
						node->outEdges[sym]->inEdges.push_back(std::make_pair(sym, node));
					}
					node = node->outEdges[sym];
				}
				else
				{
					node->outEdges[sym] = newRoot;
					newRoot->inEdges.push_back(std::make_pair(sym, node));
				}
			}
		}
		currentRoot = newRoot;
	}
	return root;
}


void Xaligner::topSort(GraphNode* root, std::vector<GraphNode*>& out)
{
	std::unordered_set<GraphNode*> visited;
	std::function<void(GraphNode*)> dfs = [&visited, &out, &dfs] (GraphNode* node)
	{
		visited.insert(node);
		for (auto nextNode : node->outEdges)
		{
			if (visited.find(nextNode.second) == visited.end())
			{
				dfs(nextNode.second);
			}
		}
		out.push_back(node);
	};
	dfs(root);
	std::reverse(out.begin(), out.end());
}

void Xaligner::alignOnGraph(const std::string& sequence, GraphNode* root, 
							int threshold, std::vector<AlignInfo>& result)
{
	//const int INS_PEN = -1;
	//const int DEL_PEN = -1;
	
	std::vector<GraphNode*> graphNodes;
	topSort(root, graphNodes);
	std::list<std::pair<GraphNode*, int>> candidates;

	std::unordered_map<GraphNode*, std::vector<int>> D;
	std::unordered_map<GraphNode*, std::vector<std::pair<GraphNode*, int>>> prev;
	for (GraphNode* node : graphNodes)
	{
		D[node].assign(sequence.length() + 1, 0);
		prev[node].assign(sequence.length() + 1, std::make_pair(nullptr, 0));
	}

	for (size_t i = 1; i < graphNodes.size(); ++i)
	{
		GraphNode* node = graphNodes[i];
		for (int n = 1; n < sequence.length() + 1; ++n)
		{
			int matchScore = std::numeric_limits<int>::min();
			GraphNode* matchNode = nullptr;
			for (auto edge : node->inEdges)
			{
				GraphNode* prevNode = edge.second;
				char sym = edge.first;
				int score = D[prevNode][n - 1] + this->subsCost(sequence[n - 1], sym);
				if (score > matchScore)
				{
					matchScore = score;
					matchNode = prevNode;
				}
			}

			int delScore = std::numeric_limits<int>::min();
			GraphNode* delNode = nullptr;
			for (auto edge : node->inEdges)
			{
				GraphNode* prevNode = edge.second;
				int score = D[prevNode][n] + _delPen;
				if (score > delScore)
				{
					delScore = score;
					delNode = prevNode;
				}
			}

			int insScore = D[node][n - 1] + _insPen;

			int score = std::max(std::max(0, matchScore), std::max(delScore, insScore));

			if (score == matchScore)
			{
				prev[node][n] = std::make_pair(matchNode, n - 1);
			}
			else if (score == delScore)
			{
				prev[node][n] = std::make_pair(delNode, n);
			}
			else if (score == insScore)
			{
				prev[node][n] = std::make_pair(node, n - 1);
			}

			D[node][n] = score;
			if (score >= threshold)
			{
				auto p = prev[node][n];
				if (D[p.first][p.second] <= score)
				{
					candidates.remove(p);
				}
				candidates.push_back(std::make_pair(node, n));
			}
		}
	}

	//backtracking
	std::vector<AlignInfo> alignments;
	for (auto cand : candidates)
	{
		int alnScore = D[cand.first][cand.second];
		std::vector<int> path;
		auto cur = cand;
		for (;;)
		{
			if (cur.first == nullptr || D[cur.first][cur.second] <= 0) break;
			path.push_back(cur.second - 1);
			cur = prev[cur.first][cur.second];
		}
		alignments.push_back(AlignInfo(path.back(), path.front(), alnScore));
	}

	if (alignments.empty()) return;

	int maxPrev = 0;
	for (int i = 1; i < alignments.size(); ++i)
	{
		if (alignments[i].start == alignments[maxPrev].start)
		{
			if (alignments[i].score > alignments[maxPrev].score)
			{
				maxPrev = i;
			}
		}
		else
		{
			result.push_back(alignments[maxPrev]);
			maxPrev = i;
		}
	}
	result.push_back(alignments[maxPrev]);
}
