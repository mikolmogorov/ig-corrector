#ifndef XALIGN_IMPL_H
#define XALIGN_IMPL_H

#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <string>

struct AlignInfo
{
	AlignInfo(int st, int e, int sc): start(st), end(e), score(sc) {}
	int start;
	int end;
	int score;
};

class Xaligner
{
public:
	Xaligner(int matchScore = 2, int missmatchScore = -2, int insPen = -1, int delPen = -1): 
		_root(nullptr),
		_matchScore(matchScore),
		_missmatchScore(missmatchScore),
		_insPen(insPen),
		_delPen(delPen)
	{}
	~Xaligner()
	{
		this->clearGraph();
	}
	void setQuerry(const std::string& querry);
	void align(const std::string& sequence, int threshold, std::vector<AlignInfo>& result);
	
private:
	int subsCost(char a, char b)
	{
		return (a == b) ? _matchScore : _missmatchScore;
	}

	struct GraphNode
	{
		std::unordered_map<char, GraphNode*> outEdges;
		std::list<std::pair<char, GraphNode*>> inEdges;
		int id;
	};

	void clearGraph();
	GraphNode* buildGraph(const std::string& sequence);
	void topSort(GraphNode* root, std::vector<GraphNode*>& out);
	void alignOnGraph(const std::string& sequence, GraphNode* root, 
					  int threshold, std::vector<AlignInfo>& result);

	GraphNode* _root;
	int _matchScore;
	int _missmatchScore;
	int _insPen;
	int _delPen;
};

#endif
