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
	counter = 0
	def __init__(self):
		self.edges = {}
		self.in_edges = []
		self.nid = Node.counter
		Node.counter += 1

def build_graph(seqence):
	table = back_trans_table()
	root = Node()
	curRoot = root
	variant = None
	for acid in seqence:
		newRoot = Node()
		if acid == "[":
			variant = ""
			continue
		elif acid == "]":
			codons = sum([table[n] for n in variant], [])
			variant = None
		elif variant is not None:
			variant += acid
			continue
		else:
			codons = table[acid]

		for codon in codons:
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
	return 2 if a == b else -2

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
					
def loc_align(seqence, root, threshold):
	INS_PENALITY = -1
	DEL_PENALITY = -1

	graph = top_sort(root)
	candidates = []
	maxScore = 0

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
			maxScore = max(maxScore, score)
			if score >= threshold:
				try:
					p = prev[node][n]
					if D[p[0]][p[1]] <= score:
						candidates.remove(p)
				except:
					pass
				candidates.append((node, n))

	alignments = []
	for c in candidates:
		aln_score = D[c[0]][c[1]]
		path = []
		n = c
		while True:
			if n[0] == None or D[n[0]][n[1]] <= 0:
				break
			path.append(n[1] - 1)
			#print ((n[1] - 1, n[0].nid)),
			n = prev[n[0]][n[1]]
		alignments.append((path[-1], path[0], aln_score))
	
	if len(alignments) == 0:
		return []
	
	result = []
	maxPrev = 0
	for i in xrange(1, len(alignments)):
		if alignments[i][0] == alignments[maxPrev][0]:
			if alignments[i][2] > alignments[maxPrev][2]:
				maxPrev = i
		else:
			result.append(alignments[maxPrev])
			maxPrev = i
	result.append(alignments[maxPrev])
	
	return result

def main():
	graph = build_graph("YYC")
	#print loc_align("TATTATTGC ", graph, 15)
	print loc_align( "ATAGATAGACGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCCGAAGTGCAGC" + 
					"TGGTGCAGTCTGGGGGAACCTTGGTGCAGGTTGGGGGTTCTCTGACACTCTCCTGTGCAGCCTCTGGATTCA" +
					"CCTTCGGAGGTTATGACATGAGCTGGGTCCGCCAGGCTCCAGGAAAGGGGCCCGAGTGGGTCTCACGTATTA" +
					"ATATGCGTGGTGTTACCACATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCA" +
					"AGGACACGGTGTTTCTGCAAATGAACAGCCTGAAATTCGAGGACTCGGCCGTG TATTATTGC ACAGGTGGGGT" +
					"CTTCGTTAGTAGCTGGTCGGGGAGCGCCTTGGATTACTGGGGCAAAGGGACAATGGTCACCGTCTCTTCAGG" +
					"CTCGAGTGCGTCTACAAAAGGCCCGTCTGTGGCGACTATAT", graph, 15)

if __name__ == "__main__":
	main()
