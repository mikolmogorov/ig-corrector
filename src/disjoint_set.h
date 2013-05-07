#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

template <class T>
struct SetNode
{
	//SetNode(): parent(0), rank(0) {}
	SetNode(const T& val): parent(0), rank(0), data(val) {}
	SetNode* parent;
	int rank;
	T data;
};

template <class T> 
SetNode<T>* makeNode(const T& val)
{
	return new SetNode<T>(val);
}

template <class T> 
SetNode<T>* findSet(SetNode<T>* elem)
{
	if (elem->parent != 0)
	{
		elem->parent = findSet(elem->parent);
		return elem->parent;
	}
	return elem;
}

template <class T> 
void unionSet(SetNode<T>* root1, SetNode<T>* root2)
{
	if (root1->rank > root2->rank)
	{
		root2->parent = root1;
	}
	else
	{
		root1->parent = root2;
		if (root1->rank == root2->rank)
		{
			++root2->rank;
		}
	}
}

#endif
