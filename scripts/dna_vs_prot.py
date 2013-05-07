#!/usr/bin/env python

from Bio.Data import CodonTable
import sys

def back_trans_table():
	table = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
	back_table = {} 
	for key in sorted(table): 
		if not table[key] in back_table:
			back_table[table[key]] = []
		back_table[table[key]].append(key)
	back_table["*"] = ["TAG", "TAA", "TGA"] 
	return back_table

class Node:
	def __init__(self):
		self.edges = {}
		self.in_edges = []

def build_graph(seqence):
	table = back_trans_table()
	root = Node()
	curRoot = root
	for acid in seqence:
		newRoot = Node()
		for codon in table[acid]:
			node = curRoot
			for i, sym in enumerate(codon):
				if i < 2:
					if not sym in node.edges:
						node.edges[sym] = Node()
						node.edges[sym].in_edges.append((sym, node))
					node = node.edges[sym]
				else:
					node.edges[sym] = newRoot
					newRoot.in_edges.append((sym, node))
		curRoot = newRoot
	return root

def subst_cost(a, b):
	return 0 if a == b else 1

def dfs(visited, path, node):
	visited.add(node)
	for next_node in node.edges.values():
		if next_node not in visited:
			dfs(visited, path, next_node)
	path.append(node)

def top_sort(root):
	visited = set()
	path = []
	dfs(visited, path, root)
	return path[::-1]
					
def compute_dist(seqence, root):
	graph = top_sort(root)
	D = {graph[0] : [i for i in xrange(len(seqence) + 1)]}

	for node in graph[1:]:
		if node not in D:
			D[node] = [0] * (len(seqence) + 1)
			
		maxInit = sys.maxint
		for _, prevNode in node.in_edges:
			maxInit = min(D[prevNode][0] + 1, maxInit)
		D[node][0] = maxInit

		for n in xrange(1, len(seqence) + 1):
			matchScore = sys.maxint
			for sym, prevNode in node.in_edges:
				score = D[prevNode][n - 1] + subst_cost(seqence[n - 1], sym)
				matchScore = min(score, matchScore)

			delScore = sys.maxint
			for _, prevNode in node.in_edges:
				score = D[prevNode][n] + 1
				delScore = min(score, delScore)

			insScore = D[node][n - 1] + 1
			D[node][n] = min(matchScore, delScore, insScore)

	print [(i,x.edges.keys()) for i, x in enumerate(graph)]
	for node in graph:
		for n in xrange(0, len(seqence) + 1):
			print D[node][n],
		print ""
	return D[graph[-1]][len(seqence)]
	
graph = build_graph("VL")
print compute_dist("TATTAC", graph)
