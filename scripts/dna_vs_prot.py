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
		#print acid
		newRoot = Node()
		for codon in table[acid]:
			node = curRoot
			for i, sym in enumerate(codon):
				#print "\t", sym
				if i < 2:
					if not sym in node.edges:
						node.edges[sym] = Node()
						node.edges[sym].in_edges.append(node)
					node = node.edges[sym]
				else:
					node.edges[sym] = newRoot
					newRoot.in_edges.append(node)
		curRoot = newRoot
	return root

def build_graph2(seqence):
	table = back_trans_table()
	root = Node()
	curRoot = root
	graph = [root]
	for acid in seqence:
		#print acid
		newRoot = Node()
		for codon in table[acid]:
			node = curRoot
			for i, sym in enumerate(codon):
				#print "\t", sym
				if i < 2:
					if not sym in node.edges:
						node.edges[sym] = Node()
						graph.append(node.edges[sym])
						node.edges[sym].in_edges.append(node)
					node = node.edges[sym]
				else:
					node.edges[sym] = newRoot
					graph.append(node.edges[sym])
					newRoot.in_edges.append(node)
		curRoot = newRoot
	return graph

def subst_cost(a, b):
	return 0 if a == b else 1

def compute_dist(seqence, root):
	D = {root : [i for i in xrange(len(seqence) + 1)]}

	counter = 2
	node_enum = {graph : 1}

	queue = [root]
	visited = []

	while len(queue) > 0:
		node = queue.pop(0)

		for sym, ch_node in node.edges.iteritems():
			if not ch_node in node_enum:
				node_enum[ch_node] = counter
				counter += 1					

			if ch_node not in visited:
				queue.append(ch_node)	
				visited.append(ch_node)
			else:
				continue

			print node_enum[node], node_enum[ch_node], sym, len(node.edges)

			assert ch_node not in D
			D[ch_node] = [0] * (len(seqence) + 1)
			
			maxInit = sys.maxint
			for prevNode in ch_node.in_edges:
				maxInit = min(D[prevNode][0] + 1, maxInit)
			D[ch_node][0] = maxInit

			for n in xrange(1, len(seqence) + 1):
				#print sym, seqence[n - 1],

				matchScore = sys.maxint
				for prevNode in ch_node.in_edges:
					#print prevNode.edges.keys(),
					symbol = [x for x in prevNode.edges if prevNode.edges[x] == ch_node][0]
					score = D[prevNode][n - 1] + subst_cost(seqence[n - 1], symbol)
					#print score,
					matchScore = min(score, matchScore)
				#print ""

				delScore = sys.maxint
				for prevNode in ch_node.in_edges:
					score = D[prevNode][n] + 1
					delScore = min(score, delScore)

				insScore = D[ch_node][n - 1] + 1

				#print matchScore, delScore, insScore, "\n"
				D[ch_node][n] = min(matchScore, delScore, insScore)
	return D[node][len(seqence)]
					
def compute_dist2(seqence, graph):
	D = {graph[0] : [i for i in xrange(len(seqence) + 1)]}

	for node in graph[1:]:
		if node not in D:
			D[node] = [0] * (len(seqence) + 1)
			
		maxInit = sys.maxint
		for prevNode in node.in_edges:
			maxInit = min(D[prevNode][0] + 1, maxInit)
		D[node][0] = maxInit

		for n in xrange(1, len(seqence) + 1):
			matchScore = sys.maxint
			for prevNode in node.in_edges:
				sym = [x for x in prevNode.edges if prevNode.edges[x] == node][0]
				score = D[prevNode][n - 1] + subst_cost(seqence[n - 1], sym)
				#print seqence[n - 1], sym
				matchScore = min(score, matchScore)

			delScore = sys.maxint
			for prevNode in node.in_edges:
				score = D[prevNode][n] + 1
				delScore = min(score, delScore)

			insScore = D[node][n - 1] + 1

			D[node][n] = min(matchScore, delScore, insScore)

	return D[graph[-1]][len(seqence)]

def print_graph(graph):
	queue = [graph]
	visited = []

	counter = 2
	node_enum = {graph : 1}
	while len(queue) > 0:
		node = queue.pop(0)
		for sym, ch_node in node.edges.iteritems():
			if not ch_node in node_enum:
				node_enum[ch_node] = counter
				counter += 1
			if ch_node not in visited:
				queue.append(ch_node)	
			visited.append(ch_node)
			print node_enum[node], node_enum[ch_node], sym
	
graph = build_graph2("YY")
print compute_dist2("TATTAT", graph)
