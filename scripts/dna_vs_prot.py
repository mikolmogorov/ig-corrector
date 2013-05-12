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
	return 2 if a == b else -1

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
					
def compute_dist(seqence, root, threshold):
	INS_PENALITY = -1
	DEL_PENALITY = -1

	graph = top_sort(root)
	candidates = []

	D = {node : [0] * (len(seqence) + 1) for node in graph}
	prev = {node : [(None, 0)] * (len(seqence) + 1) for node in graph}

	for node in graph[1:]:
		for n in xrange(1, len(seqence) + 1):
			matchScore = -sys.maxint
			matchNode = None
			for sym, prevNode in node.in_edges:
				score = D[prevNode][n - 1] + subst_cost(seqence[n - 1], sym)
				if score > matchScore:
					matchScore = score
					matchNode = prevNode

			delScore = -sys.maxint
			delNode = None
			for _, prevNode in node.in_edges:
				score = D[prevNode][n] + DEL_PENALITY
				if score > delScore:
					delScore = score
					delNode = prevNode

			insScore = D[node][n - 1] + INS_PENALITY

			score = max(0, matchScore, delScore, insScore)

			if score == matchScore:
				prev[node][n] = (matchNode, n - 1)
			elif score == delScore:
				prev[node][n] = (delNode, n)
			elif score == insScore:
				prev[node][n] = (node, n - 1)

			D[node][n] = score
			if score >= threshold:
				try:
					candidates.remove(prev[node][n])
				except:
					pass
				candidates.append((node, n))

	alignments = []
	for c in candidates:
		path = []
		n = c
		while True:
			if n[0] == None or D[n[0]][n[1]] <= 0:
				break
			path.append(n[1] - 1)
			n = prev[n[0]][n[1]]
		alignments.append((path[-1], path[0]))

	########
	#print [(i,x.edges.keys()) for i, x in enumerate(graph)]
	#for node in graph:
	#	for n in xrange(0, len(seqence) + 1):
	#		print "%2d" % D[node][n],
	#	print ""
	#######
	return alignments
	
graph = build_graph("YYC")
print compute_dist( "ATAGATAGACGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCCGAAGTGCAGC" + 
					"TGGTGCAGTCTGGGGGAACCTTGGTGCAGGTTGGGGGTTCTCTGACACTCTCCTGTGCAGCCTCTGGATTCA" +
					"CCTTCGGAGGTTATGACATGAGCTGGGTCCGCCAGGCTCCAGGAAAGGGGCCCGAGTGGGTCTCACGTATTA" +
					"ATATGCGTGGTGTTACCACATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCA" +
					"AGGACACGGTGTTTCTGCAAATGAACAGCCTGAAATTCGAGGACTCGGCCGTGTATACTGTACAGGTGGGGT" +
					"CTTCGTTAGTAGCTGGTCGGGGAGCGCCTTGGATTACTGGGGCAAAGGGACAATGGTCACCGTCTCTTCAGG" +
					"CTCGAGTGCGTCTACAAAAGGCCCGTCTGTGGCGACTATAT", graph, 16)
