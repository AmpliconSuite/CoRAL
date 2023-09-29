import numpy as np
from breakpoint_graph import *
from cycle_decomposition import *

if __name__ == '__main__':
	g = BreakpointGraph()
	fp = open("/ribosome/projects/Nanopore/LR_output/hybrid_new/BT474_amplicon1_cycles.txt", 'r')
	#fp = open("/ribosome/projects/Nanopore/LR_output/hybrid_new/COLO320DM_Wu2019_amplicon4_cycles.txt", 'r')
	#fp = open("/ribosome/projects/Nanopore/LR_output/hybrid_new/SKBR3_550_amplicon2_cycles.txt", 'r')
	
	for line in fp:
		s = line.strip().split()
		if s[0] == 'Interval':
			g.amplicon_intervals.append([s[2], int(s[3]), int(s[4])])
			g.add_endnode((s[2], int(s[3]), '-'))
			g.add_endnode((s[2], int(s[4]), '+'))
	fp.close()


	fp = open("/ribosome/projects/Nanopore/LR_output/hybrid_new/BT474_amplicon1_graph.txt", 'r')
	#fp = open("/ribosome/projects/Nanopore/LR_output/hybrid_new/COLO320DM_Wu2019_amplicon4_graph.txt", 'r')
	#fp = open("/ribosome/projects/Nanopore/LR_output/hybrid_new/SKBR3_550_amplicon2_graph.txt", 'r')
	max_cn = 0.0
	for line in fp:
		s = line.strip().split("\t")
		if s[0] == 'sequence':
			chr = s[1][: s[1].find(':')]
			l = int(s[1][s[1].find(':') + 1: -1])
			r = int(s[2][s[2].find(':') + 1: -1])
			g.add_node((chr, l, '-'))
			g.add_node((chr, r, '+'))
			g.add_sequence_edge(chr, l, r, int(s[4]), 'd', int(s[5]), 0, float(s[3]))
			if float(s[3]) > max_cn:
				max_cn = float(s[3])		
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
	g.max_cn = max_cn
	g.del_discordant_endnodes()
	print (len(g.sequence_edges), len(g.concordant_edges), len(g.discordant_edges), len(g.source_edges))
	node_order = dict()
	ni = 0
	for node in g.nodes.keys():
		node_order[node] = ni
		ni += 1
	s, f, r = max_flow(g, node_order, 8, "BT474_amplicon1")
	
	f.merge_edges()
	print (len(f.sequence_edges), len(f.concordant_edges), len(f.discordant_edges), len(f.source_edges))
	
	clean_up_residual_network(r, 0.1)
	r.merge_edges()
	print (len(r.sequence_edges), len(r.concordant_edges), len(r.discordant_edges), len(r.source_edges))
	
	node_order_f = dict()
	ni = 0
	for node in f.nodes.keys():
		node_order_f[node] = ni
		ni += 1
	
	node_order_r = dict()
	ni = 0
	for node in r.nodes.keys():
		node_order_r[node] = ni
		ni += 1
	
	#s1 = minimize_cycles_max_flow(f, node_order_f, 2, 8, "BT474_amplicon1")
	#total_weights_f = sum([seq[-1] * seq[-2] for seq in f.sequence_edges])
	total_weights_r = sum([seq[-1] * seq[-2] for seq in r.sequence_edges])
	avg_cn = np.average([seq[-1] for seq in r.sequence_edges], weights = [seq[-2] for seq in r.sequence_edges])
	max_cn = max([seq[-1] for seq in r.sequence_edges if seq[-2] > 50000])
	print (avg_cn, max_cn)
	"""
	for seq in r.sequence_edges:
		if seq[-1] == 0.0:
			print ('1', seq)
		else:
			print ('2', seq)
	"""
	#print ('w=', total_weights_f)
	print ('w=', total_weights_r)
	
	#maximize_weights_greedy(f, total_weights_f, node_order_f, [], 0.01, 2, 0.9, 0.1, 8, "BT474_amplicon1_f")
	maximize_weights_greedy_(r, total_weights_r, node_order_r, [], 0.01, 4, 0.9, 0.1, 16, "BT474_amplicon1_r")
	#print (f.sequence_edges)
	#print (f.concordant_edges)
	#print (f.discordant_edges)



