from breakpoint_graph import *
from cycle_decomposition import *

if __name__ == '__main__':
	g = BreakpointGraph()
	fp = open("/ribosome/projects/Nanopore/LR_output/hybrid_new/BT474_amplicon1_cycles.txt", 'r')
	for line in fp:
		s = line.strip().split()
		if s[0] == 'Interval':
			g.amplicon_intervals.append([s[2], int(s[3]), int(s[4])])
			g.add_endnode((s[2], int(s[3]), '-'))
			g.add_endnode((s[2], int(s[4]), '+'))
	fp.close()


	fp = open("/ribosome/projects/Nanopore/LR_output/hybrid_new/BT474_amplicon1_graph.txt", 'r')
	for line in fp:
		s = line.strip().split("\t")
		if s[0] == 'sequence':
			chr = s[1][: s[1].find(':')]
			l = int(s[1][s[1].find(':') + 1: -1])
			r = int(s[2][s[2].find(':') + 1: -1])
			g.add_node((chr, l, '-'))
			g.add_node((chr, r, '+'))
			g.add_sequence_edge(chr, l, r, int(s[4]), 'd', int(s[5]), 0, float(s[3]))		
		elif s[0] == 'concordant':
			t = s[1].split('->')
			t0_ = t[0].split(':')
			t1_ = t[1].split(':')
			assert t0_[1][-1] == '+' and t1_[1][-1] == '-' and t0_[0] == t1_[0] and int(t1_[1][:-1]) == int(t0_[1][:-1]) + 1
			g.add_concordant_edge(t0_[0], int(t0_[1][:-1]), t0_[1][-1], t1_[0], int(t1_[1][:-1]), t1_[1][-1], 
						int(s[3]), 'd', int(s[4]), set([]), float(s[2]))
		elif s[0] == 'discordant':
			t = s[1].split('->')
			t0_ = t[0].split(':')
			t1_ = t[1].split(':')
			g.add_discordant_edge(t0_[0], int(t0_[1][:-1]), t0_[1][-1], t1_[0], int(t1_[1][:-1]), t1_[1][-1], 
						max(int(s[3]), 0), 'd', 0.0, int(s[4]), set([]), float(s[2]))
		elif s[0] == 'source':
			t = s[1].split('->')
			t0_ = t[0].split(':')
			t1_ = t[1].split(':')
			g.add_source_edge(t1_[0], int(t1_[1][:-1]), t1_[1][-1], int(round(float(s[3]))), 'd', 0.0, 0.0, float(s[2]))

	fp.close()
	print (len(g.endnodes))
	g.del_discordant_endnodes()
	print (len(g.endnodes))
	node_order = dict()
	ni = 0
	for node in g.nodes.keys():
		node_order[node] = ni
		ni += 1
	s, f, r = max_flow(g, node_order, 8, "BT474_amplicon1")
	print (f.sequence_edges)
	print (f.concordant_edges)
	print (f.discordant_edges)



