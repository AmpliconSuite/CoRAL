"""
Functions used for cycle decomposition
"""
import math
import copy
import gurobipy as gp
from gurobipy import GRB

from breakpoint_graph import *


def max_flow(g, node_order, num_threads = 16, model_prefix = ""):
	lseg = len(g.sequence_edges)
	lc = len(g.concordant_edges)
	ld = len(g.discordant_edges)
	lsrc = len(g.source_edges)
	nedges = lseg + lc + ld + lsrc + len(g.endnodes)
	endnode_list = [node for node in g.endnodes.keys()]

	m = gp.Model(model_prefix + "_max_flow")

	# y: flow values (directed)
	y1_names = []
	y2_names = []
	for ei in range(nedges):
		y1_names.append("y1_" + str(ei))
		y2_names.append("y2_" + str(ei))
	y1 = m.addVars(nedges, vtype = GRB.CONTINUOUS, lb = 0.0, name = y1_names) #small -> large
	y2 = m.addVars(nedges, vtype = GRB.CONTINUOUS, lb = 0.0, name = y2_names) #large -> small

	# Objective: maximize the flow leaving s
	obj = gp.LinExpr(0.0)
	for yi in range(lsrc + len(g.endnodes)):
		obj += (y1[lseg + lc + ld + yi])
	m.setObjective(obj, GRB.MAXIMIZE)

	# Capacity constraints
	for seqi in range(lseg):
		m.addConstr(y1[seqi] + y2[seqi] <= g.sequence_edges[seqi][-1])
	for ci in range(lc):
		m.addConstr(y1[lseg + ci] + y2[lseg + ci] <= g.concordant_edges[ci][-1])
	for di in range(ld):
		m.addConstr(y1[lseg + lc + di] + y2[lseg + lc + di] <= g.discordant_edges[di][-1])
	for srci in range(lsrc):
		m.addConstr(y1[lseg + lc + ld + srci] + y2[lseg + lc + ld + srci] <= g.source_edges[srci][-1])
	for node in g.endnodes.keys():
		sseg = g.sequence_edges[g.nodes[node][0][0]]
		yi = endnode_list.index(node)
		m.addConstr(y1[lseg + lc + ld + lsrc + yi] + y2[lseg + lc + ld + lsrc + yi] <= sseg[-1])

	# Flow conservation constraints
	for node in g.endnodes.keys():
		seqi = g.nodes[node][0][0]
		sseg = g.sequence_edges[seqi]
		"""
		node_ = (sseg[0], sseg[1], '-')
		if node_ == node:
			node_ = (sseg[0], sseg[2], '+')
		"""
		yi = endnode_list.index(node)
		m.addConstr(y2[seqi] + y1[seqi] == y1[lseg + lc + ld + lsrc + yi] + y2[lseg + lc + ld + lsrc + yi])
		"""
		if node_order[node] < node_order[node_]:
			m.addConstr(y2[seqi] + y1[lseg + lc + ld + lsrc + yi] == y1[seqi] + y2[lseg + lc + ld + lsrc + yi])
		else:
			m.addConstr(y1[seqi] + y1[lseg + lc + ld + lsrc + yi] == y2[seqi] + y2[lseg + lc + ld + lsrc + yi])
		"""
	expr_y_endnode = gp.LinExpr(0.0)
	for yi in range(lsrc + len(g.endnodes)):
		expr_y_endnode += (y1[lseg + lc + ld + yi])
		expr_y_endnode -= (y2[lseg + lc + ld + yi])
	m.addConstr(expr_y_endnode == 0.0)
	"""
	for di in range(ld):
		de = g.discordant_edges[di]
		if de[0] == de[3] and de[1] == de[4] and de[2] == de[5]:
			node = (de[0], de[1], de[2])
			expr_y_selfloop = gp.LinExpr(0.0)
			for seqi in g.nodes[node][0]:
				expr_y_selfloop += (y1[seqi] + y2[seqi])
			for ci in g.nodes[node][1]:
				expr_y_selfloop -= (y1[lseg + ci] + y2[lseg + ci])
			for di_ in g.nodes[node][2]:
				if di_ == di:
					expr_y_selfloop -= 2 * (y1[lseg + lc + di] + y2[lseg + lc + di])
				else:
					expr_y_selfloop -= (y1[lseg + lc + di] + y2[lseg + lc + di])
			m.addConstr(expr_y_selfloop == 0.0)
	"""
	for node in g.nodes.keys():
		if node not in g.endnodes:
			expr_y = gp.LinExpr(0.0)
			for seqi in g.nodes[node][0]:
				sseg = g.sequence_edges[seqi]
				node_ = (sseg[0], sseg[1], '-')
				if node_ == node:
					node_ = (sseg[0], sseg[2], '+')
				if node_order[node] < node_order[node_]:
					expr_y += y2[seqi]
					expr_y += y1[seqi]
				else:
					expr_y += y1[seqi]
					expr_y += y2[seqi]
			for ci in g.nodes[node][1]:
				ce = g.concordant_edges[ci]
				node_ = (ce[0], ce[1], ce[2])
				if node_ == node:
					node_ = (ce[3], ce[4], ce[5])
				if node_order[node] < node_order[node_]:
					expr_y -= y2[lseg + ci]
					expr_y -= y1[lseg + ci]
				else:
					expr_y -= y1[lseg + ci]
					expr_y -= y2[lseg + ci]
			for di in g.nodes[node][2]:
				de = g.discordant_edges[di]
				node_ = (de[0], de[1], de[2])
				if node_ == node:
					node_ = (de[3], de[4], de[5])
				if node_order[node] < node_order[node_]:
					expr_y -= y2[lseg + lc + di]
					expr_y -= y1[lseg + lc + di]
				elif node_order[node] > node_order[node_]:
					expr_y -= y1[lseg + lc + di]
					expr_y -= y2[lseg + lc + di]
				else:
					expr_y -= 2 * y1[lseg + lc + di]
					expr_y -= 2 * y2[lseg + lc + di]
			for srci in g.nodes[node][3]:
				expr_y -= y1[lseg + lc + ld + srci]
				expr_y -= y2[lseg + lc + ld + srci]
			m.addConstr(expr_y == 0.0)
			
	m.setParam(GRB.Param.Threads, num_threads)
	m.write(model_prefix + "_max_flow.lp") 
	m.optimize()
	flow_network = BreakpointGraph()
	flow_network.max_cn = g.max_cn
	residual_network = copy.deepcopy(g)
	if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
		print('Optimization was stopped with status %d' % m.Status)
		return GRB.INFEASIBLE, flow_network, residual_network
	else:
		sol_y1 = m.getAttr('X', y1)
		sol_y2 = m.getAttr('X', y2)
		for ni in range(len(endnode_list)):
			node = endnode_list[ni]
			#sseg = g.sequence_edges[g.nodes[node][0][0]]
			if sol_y1[lseg + lc + ld + lsrc + ni] + sol_y2[lseg + lc + ld + lsrc + ni] > 0:
				flow_network.add_endnode(node)
				#if node == ('chr1', 203349098, '-'):
				#	print (sol_y1[lseg + lc + ld + lsrc + ni], sol_y2[lseg + lc + ld + lsrc + ni])
		for yi in range(len(sol_y1)):
			if sol_y1[yi] + sol_y2[yi] > 0:
				if yi < lseg:
					sseg = g.sequence_edges[yi]
					residual_network.sequence_edges[yi][-1] -= (sol_y1[yi] + sol_y2[yi])
					flow_network.add_node((sseg[0], sseg[1], '-'))
					flow_network.add_node((sseg[0], sseg[2], '+'))
					flow_network.add_sequence_edge(sseg[0], sseg[1], sseg[2], cn = (sol_y1[yi] + sol_y2[yi]))
				elif yi < lseg + lc:
					ce = g.concordant_edges[yi - lseg]
					residual_network.concordant_edges[yi - lseg][-1] -= (sol_y1[yi] + sol_y2[yi])
					flow_network.add_concordant_edge(ce[0], ce[1], ce[2], ce[3], ce[4], ce[5], cn = (sol_y1[yi] + sol_y2[yi]))
				elif yi < lseg + lc + ld:
					de = g.discordant_edges[yi - lseg - lc]
					residual_network.discordant_edges[yi - lseg - lc][-1] -= (sol_y1[yi] + sol_y2[yi])
					#print (de, sol_y1[yi], sol_y2[yi])
					flow_network.add_discordant_edge(de[0], de[1], de[2], de[3], de[4], de[5], cn = (sol_y1[yi] + sol_y2[yi]))
				elif yi < lseg + lc + ld + lsrc:
					srce = g.source_edges[yi - lseg - lc - ld]
					residual_network.source_edges[yi - lseg - lc - ld][-1] -= (sol_y1[yi] + sol_y2[yi])
					flow_network.add_source_edge(srce[3], srce[4], srce[5], cn = (sol_y1[yi] + sol_y2[yi]))
		#for node in endnode_list:
		#	yi = endnode_list.index(node)
		#	print (node, sol_y1[lseg + lc + ld + lsrc + yi], sol_y2[lseg + lc + ld + lsrc + yi])
		#for node in flow_network.endnodes.keys():
		#	if node not in flow_network.nodes:
		#		print ('N', node)
		return m.Status, flow_network, residual_network		


def clean_up_residual_network(r, resolution = 0.1):
	lseg = len(r.sequence_edges)
	lc = len(r.concordant_edges)
	ld = len(r.discordant_edges)
	lsrc = len(r.source_edges)
	for seqi in range(lseg)[::-1]:
		sseg = r.sequence_edges[seqi]
		if sseg[-1] < resolution:
			node1 = (sseg[0], sseg[1], '-')
			node2 = (sseg[0], sseg[2], '+')
			del r.sequence_edges[seqi]
			if node1 in r.nodes:
				del r.nodes[node1]
			if node2 in r.nodes:
				del r.nodes[node2]
			if node1 in r.endnodes:
				del r.endnodes[node1]
			if node2 in r.endnodes:
				del r.endnodes[node2]
	for seqi in range(len(r.sequence_edges)):
		sseg = r.sequence_edges[seqi]
		node1 = (sseg[0], sseg[1], '-')
		node2 = (sseg[0], sseg[2], '+')
		r.nodes[node1][0] = [seqi]
		r.nodes[node2][0] = [seqi]
		r.nodes[node1][1] = []
		r.nodes[node2][1] = []
		r.nodes[node1][2] = []
		r.nodes[node2][2] = []
		r.nodes[node1][3] = []
		r.nodes[node2][3] = []
	for ci in range(lc)[::-1]:
		ce = r.concordant_edges[ci]
		node1 = (ce[0], ce[1], ce[2])
		node2 = (ce[3], ce[4], ce[5])
		if ce[-1] < resolution or node1 not in r.nodes or node2 not in r.nodes:
			del r.concordant_edges[ci]
	for ci in range(len(r.concordant_edges)):
		ce = r.concordant_edges[ci]
		node1 = (ce[0], ce[1], ce[2])
		node2 = (ce[3], ce[4], ce[5])
		r.nodes[node1][1] = [ci]
		r.nodes[node2][1] = [ci]
	for di in range(ld)[::-1]:
		de = r.discordant_edges[di]
		node1 = (de[0], de[1], de[2])
		node2 = (de[3], de[4], de[5])
		if de[-1] < resolution or node1 not in r.nodes or node2 not in r.nodes:
			del r.discordant_edges[di]
	for di in range(len(r.discordant_edges)):
		de = r.discordant_edges[di]
		node1 = (de[0], de[1], de[2])
		node2 = (de[3], de[4], de[5])
		r.nodes[node1][2].append(di)
		r.nodes[node2][2].append(di)
	for srci in range(lsrc)[::-1]:
		srce = r.source_edges[srci]
		node1 = (srce[3], srce[4], srce[5])
		if srce[-1] < resolution or node1 not in r.nodes or node2 not in r.nodes:
			del r.source_edges[srci]
	for srci in range(len(r.source_edges)):
		srce = r.source_edges[srci]
		node1 = (srce[3], srce[4], srce[5])
		r.nodes[node1][3].append(srci)
	"""
	for node in r.nodes.keys():
		if r.nodes[node][0][0] > len(r.sequence_edges):
			print ('s', node, r.nodes[node])
		if r.nodes[node][1][0] > len(r.concordant_edges):
			print ('c', node, r.nodes[node])
	"""


def minimize_cycles_max_flow(g, node_order, max_seq_repeat = 2, num_threads = 16, model_prefix = ""):
	"""
	Given a breakpoint graph corresponding the max flow network, 
	compute a minimum number of s-t walks explaining all the CN weights
	"""
	lseg = len(g.sequence_edges)
	lc = len(g.concordant_edges)
	ld = len(g.discordant_edges)
	lsrc = len(g.source_edges)
	nedges = lseg + lc + ld + 2 * lsrc + 2 * len(g.endnodes)
	nnodes = len(g.nodes)
	endnode_list = [node for node in g.endnodes.keys()]

	m = gp.Model(model_prefix + "cycle_decomposition_max_flow")

	k = lsrc + len(endnode_list) + ld

	# z[i]: indicating whether cycle or path i exists
	z = m.addVars(k, vtype = GRB.BINARY, name = ["z" + str(i) for i in range(k)])
		
	# w[i]: the weight of cycle or path i, continuous variable
	w = m.addVars(k, lb = 0.0, ub = g.max_cn, vtype = GRB.CONTINUOUS, name = ["w" + str(i) for i in range(k)])
		
	# Relationship between w[i] and z[i]
	for i in range(k):
		m.addConstr(w[i] <= z[i] * g.max_cn)

	# x: the number of times an edge occur in cycle or path i
	x_names = []
	for ei in range(nedges):
		for i in range(k):
			x_names.append("x" + str(ei) + "," + str(i))
	x = m.addVars(k * nedges, lb = 0.0, ub = 10.0, vtype = GRB.INTEGER, name = x_names)
	
	# Objective: minimize the total number of cycles
	obj = gp.LinExpr(0.0)
	for i in range(k):
		obj += z[i]
	m.setObjective(obj, GRB.MINIMIZE)

	# Eulerian constraint
	for node in g.nodes.keys():
		if node in endnode_list:
			for i in range(k):
				m.addConstr(x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + i] + \
						x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i] \
						== x[g.nodes[node][0][0] * k + i])
		else:
			for i in range(k):
				ec_expr = gp.LinExpr(0.0)
				for seqi in g.nodes[node][0]:
					ec_expr += x[seqi * k + i]
				for ci in g.nodes[node][1]:
					ec_expr -= x[(lseg + ci) * k + i]
				for di in g.nodes[node][2]:
					de = g.discordant_edges[di]
					if de[0] == de[3] and de[1] == de[4] and de[2] == de[5]:
						ec_expr -= 2 * x[(lseg + lc + di) * k + i]
					else:
						ec_expr -= x[(lseg + lc + di) * k + i]
				for srci in g.nodes[node][3]:
					ec_expr -= x[(lseg + lc + ld + 2 * srci) * k + i] # connected to s
					ec_expr -= x[(lseg + lc + ld + 2 * srci) * k + k + i] # connected to t
				m.addConstr(ec_expr == 0.0)
	for i in range(k):
		path_expr = gp.LinExpr(0.0)
		for ni in range(len(endnode_list)):
			path_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + i] # (s, v)
			path_expr -= x[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + k + i] # (v, t)
		for srci in range(lsrc): 
			path_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
			path_expr -= x[(lseg + lc + ld + 2 * srci) * k + k + i] # (v, t)
		m.addConstr(path_expr == 0.0)

	# CN constraint
	for seqi in range(lseg):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[seqi * k + i]
		m.addQConstr(cn_expr == g.sequence_edges[seqi][-1])
	for ci in range(lc):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[(lseg + ci) * k + i]
		m.addQConstr(cn_expr == g.concordant_edges[ci][-1])
	for di in range(ld):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[(lseg + lc + di) * k + i]
		m.addQConstr(cn_expr == g.discordant_edges[di][-1])
	for srci in range(lsrc):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + i]
			cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + k + i]
		m.addQConstr(cn_expr == g.source_edges[srci][-1])
			
	# Occurrence of breakpoints in each cycle/path
	for i in range(k):
		for seqi in range(lseg):
			m.addConstr(x[seqi * k + i] <= max_seq_repeat)
			
	# Decompose into s-t paths/walks
	for i in range(k):
		cycle_expr = gp.LinExpr(0.0)
		for ni in range(len(endnode_list)):
			cycle_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + i] # (s, v)
		for srci in range(lsrc): 
			cycle_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
		m.addConstr(cycle_expr <= 1.0)
		
	# d: BFS/spanning tree order of the nodes in decomposition i
	d_names = []
	for ni in range(nnodes):
		for i in range(k):
			d_names.append("d" + str(ni) + "," + str(i))
	for i in range(k):
		d_names.append("ds," + str(i))
	for i in range(k):
		d_names.append("dt," + str(i))
	d = m.addVars(k * (nnodes + 2), lb = 0.0, ub = nnodes + 2, vtype = GRB.INTEGER, name = d_names) # including s, t, at the end
			
	# y: spanning tree indicator (directed)
	y1_names = []
	y2_names = []
	for ei in range(nedges):
		for i in range(k):
			y1_names.append("y1-" + str(ei) + "," + str(i))
			y2_names.append("y2-" + str(ei) + "," + str(i))
	y1 = m.addVars(k * nedges, vtype = GRB.BINARY, name = y1_names) #small -> large
	y2 = m.addVars(k * nedges, vtype = GRB.BINARY, name = y2_names) #large -> small
			
	# Relationship between y and z:
	for i in range(k):
		for j in range(nedges):
			m.addConstr(y1[j * k + i] <= z[i])
			m.addConstr(y2[j * k + i] <= z[i])

	# Relationship between x, y and d
	for i in range(k * nedges):
		m.addConstr(y1[i] + y2[i] <= x[i])
	for i in range(k):
		for di in range(ld):
			de = g.discordant_edges[di]
			if de[0] == de[3] and de[1] == de[4] and de[2] == de[5]: # exclude self loops
				m.addConstr(y1[(lseg + lc + di) * k + i] == 0)
				m.addConstr(y2[(lseg + lc + di) * k + i] == 0)
	for i in range(k):
		t_expr_x = gp.LinExpr(0.0)
		t_expr_y = gp.LinExpr(0.0)
		t_expr_yd = gp.QuadExpr(0.0)
		for ni in range(len(endnode_list)):
			node = endnode_list[ni]
			t_expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + k + i]
			t_expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + k + i]
			t_expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + k + i] * \
							(d[k * nnodes + k + i] - d[k * node_order[node] + i]) # node -> t
			expr_x = gp.LinExpr(0.0)
			expr_y = gp.LinExpr(0.0)
			expr_yd = gp.QuadExpr(0.0)
			#print (node, g.nodes)
			for seqi in g.nodes[node][0]:
				sseg = g.sequence_edges[seqi]
				node_ = (sseg[0], sseg[1], '-')
				if node_ == node:
					node_ = (sseg[0], sseg[2], '+')
				expr_x += x[seqi * k + i]
				if node_order[node_] <= node_order[node]:
					expr_y += y1[seqi * k + i]
					expr_yd += y1[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				else:
					expr_y += y2[seqi * k + i]
					expr_yd += y2[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						
				expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + i] # from s
				expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + k + i] # to t
				expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + i] # from s
				expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * ni) * k + i] * \
								(d[k * node_order[node] + i] - d[k * nnodes + i])
			m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
			m.addConstr(expr_y <= 1.0)
			m.addConstr(expr_y * nedges * k >= expr_x)
			m.addConstr(expr_yd * nedges * k >= expr_x)
			
		for srci in range(lsrc):
			srce = g.source_edges[srci]
			t_expr_x += x[(lseg + lc + ld + 2 * srci) * k + k + i]
			t_expr_y += y1[(lseg + lc + ld + 2 * srci) * k + k + i]
			t_expr_yd += y1[(lseg + lc + ld + 2 * srci) * k + k + i] * \
					(d[k * nnodes + k + i] - d[k * node_order[(srce[3], srce[4], srce[5])] + i])
		m.addConstr(t_expr_x * (nnodes + 2) >= d[k * nnodes + k + i])
		m.addConstr(t_expr_y <= 1.0)
		m.addConstr(t_expr_y * nedges * k >= t_expr_x)
		m.addConstr(t_expr_yd >= t_expr_x)
			
		for node in g.nodes.keys():
			if node not in endnode_list:
				expr_x = gp.LinExpr(0.0)
				expr_y = gp.LinExpr(0.0)
				expr_yd = gp.QuadExpr(0.0)
				for seqi in g.nodes[node][0]:
					sseg = g.sequence_edges[seqi]
					node_ = (sseg[0], sseg[1], '-')
					if node_ == node:
						node_ = (sseg[0], sseg[2], '+')
					expr_x += x[seqi * k + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[seqi * k + i]
						expr_yd += y1[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[seqi * k + i]
						expr_yd += y2[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				for ci in g.nodes[node][1]:
					ce = g.concordant_edges[ci]
					node_ = (ce[0], ce[1], ce[2])
					if node_ == node:
						node_ = (ce[3], ce[4], ce[5])
					expr_x += x[(lseg + ci) * k + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[(lseg + ci) * k + i]
						expr_yd += y1[(lseg + ci) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[(lseg + ci) * k + i]
						expr_yd += y2[(lseg + ci) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				for di in g.nodes[node][2]:
					de = g.discordant_edges[di]
					node_ = (de[0], de[1], de[2])
					if node_ == node:
						node_ = (de[3], de[4], de[5])
					expr_x += x[(lseg + lc + di) * k + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[(lseg + lc + di) * k + i]
						expr_yd += y1[(lseg + lc + di) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[(lseg + lc + di) * k + i]
						expr_yd += y2[(lseg + lc + di) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				for srci in g.nodes[node][3]:
					expr_x += x[(lseg + lc + ld + 2 * srci) * k + i]
					expr_x += x[(lseg + lc + ld + 2 * srci) * k + k + i]
					expr_y += y1[(lseg + lc + ld + 2 * srci) * k + i]
					expr_yd += y1[(lseg + lc + ld + 2 * srci) * k + i] * \
							(d[k * node_order[node] + i] - d[k * nnodes + i])
				m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
				m.addConstr(expr_y <= 1.0)
				m.addConstr(expr_y * nedges * k >= expr_x)
				m.addConstr(expr_yd * nedges * k >= expr_x)
			
	m.setParam(GRB.Param.Threads, num_threads)
	m.setParam(GRB.Param.NonConvex, 2)
	m.setParam(GRB.Param.TimeLimit, max(7200, ld * 300)) # each breakpoint edge is assigned 5 minutes 
	m.write(model_prefix + "_decomp_maxflow_model.lp")
	m.optimize()
	if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
		print('Optimization was stopped with status %d' % m.Status)
		return GRB.INFEASIBLE
	else:
		return m.Status


def maximize_weights_greedy(g, total_weights, node_order, pc_list, alpha = 0.01,
			max_seq_repeat = 2, p_total_weight = 0.9, resolution = 0.1, num_threads = 16, model_prefix = ""):
	lseg = len(g.sequence_edges)
	lc = len(g.concordant_edges)
	ld = len(g.discordant_edges)
	lsrc = len(g.source_edges)
	print(lseg, lc, ld, lsrc)
	nedges = lseg + lc + ld + 2 * lsrc + 2 * len(g.endnodes)
	nnodes = len(g.nodes)
	endnode_list = [node for node in g.endnodes.keys()]

	remaining_weights = total_weights
	unsatisfied_pc = [i for i in range(len(pc_list))]
	remaining_CN = dict()
	for seqi in range(lseg):
		remaining_CN[('s', seqi)] = g.sequence_edges[seqi][-1]
	for ci in range(lc):
		remaining_CN[('c', ci)] = g.concordant_edges[ci][-1]
	for di in range(ld):
		remaining_CN[('d', di)] = g.discordant_edges[di][-1]
	for srci in range(lsrc):
		remaining_CN[('src', srci)] = g.source_edges[srci][-1]

	next_w = resolution * 1.1
	cycle_id = 0
	num_unsatisfied_pc = len(pc_list)
	while next_w >= resolution and (remaining_weights > (1.0 - p_total_weight) * total_weights or num_unsatisfied_pc > math.floor(0.1 * len(pc_list))):
		pp = 1.0
		if alpha > 0 and num_unsatisfied_pc > 0:
			pp = alpha * remaining_weights / num_unsatisfied_pc # multi - objective optimization parameter
		print ("Cycle id = ", cycle_id)
		print ("Remaining weights = ", remaining_weights, total_weights)
		print ("Num unsatisfied path constraints = ", num_unsatisfied_pc, len(pc_list))
		print ("Path constraints factor = ", pp)

		# Gurobi model
		m = gp.Model(model_prefix + "_cycle_decomposition_greedy_" + str(cycle_id + 1))

		# z[i]: indicating whether cycle or path i exists
		z = m.addVars(1, vtype = GRB.BINARY, name = ["z0"])

		# w[i]: the weight of cycle or path i, continuous variable
		w = m.addVars(1, lb = 0.0, ub = g.max_cn, vtype = GRB.CONTINUOUS, name = ["w0"])

		# Relationship between w[i] and z[i]
		m.addConstr(w[0] <= z[0] * g.max_cn)
		m.addConstr(w[0] >= z[0] * resolution)

		# x: the number of times an edge occur in cycle or path i
		x_names = []
		for ei in range(nedges):
			x_names.append("x" + str(ei))
		x = m.addVars(nedges, lb = 0.0, ub = 10.0, vtype = GRB.INTEGER, name = x_names)

		# Subpath constraints
		r_names = []
		for pi in range(len(pc_list)):
			r_names.append("r" + str(pi))
		r = m.addVars(len(pc_list), vtype = GRB.BINARY, name = r_names)
			
		# Objective: maximize the total weight + num subpath constraints satisfied
		obj = gp.QuadExpr(0.0)
		for seqi in range(lseg):
			obj += (x[seqi] * w[0] * g.sequence_edges[seqi][-2])
		for pi in range(len(pc_list)):
			if unsatisfied_pc[pi] >= 0: 
				obj += (r[pi] * max(pp, 1.0))
		m.setObjective(obj, GRB.MAXIMIZE)

		# Eulerian constraint
		for node in g.nodes.keys():
			if node in endnode_list:
				m.addConstr(x[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)] + \
						x[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node) + 1] \
						== x[g.nodes[node][0][0]])
			else:
				ec_expr = gp.LinExpr(0.0)
				for seqi in g.nodes[node][0]:
					ec_expr += x[seqi]
				for ci in g.nodes[node][1]:
					ec_expr -= x[lseg + ci]
				for di in g.nodes[node][2]:
					de = g.discordant_edges[di]
					if de[0] == de[3] and de[1] == de[4] and de[2] == de[5]:
						ec_expr -= 2 * x[lseg + lc + di]
					else:
						ec_expr -= x[lseg + lc + di]
				for srci in g.nodes[node][3]:
					ec_expr -= x[lseg + lc + ld + 2 * srci] # connected to s
					ec_expr -= x[lseg + lc + ld + 2 * srci + 1] # connected to t
				m.addConstr(ec_expr == 0.0)
		path_expr = gp.LinExpr(0.0)
		for ni in range(len(endnode_list)):
			path_expr += x[lseg + lc + ld + 2 * lsrc + 2 * ni] # (s, v)
			path_expr -= x[lseg + lc + ld + 2 * lsrc + 2 * ni + 1] # (v, t)
		for srci in range(lsrc): 
			path_expr += x[lseg + lc + ld + 2 * srci] # (s, v)
			path_expr -= x[lseg + lc + ld + 2 * srci + 1] # (v, t)
		m.addConstr(path_expr == 0.0)

		# CN constraint
		for seqi in range(lseg):
			m.addQConstr(w[0] * x[seqi] <= remaining_CN[('s', seqi)])
		for ci in range(lc):
			m.addQConstr(w[0] * x[lseg + ci] <= remaining_CN[('c', ci)])
		for di in range(ld):
			m.addQConstr(w[0] * x[lseg + lc + di] <= remaining_CN[('d', di)])
			if remaining_CN[('d', di)] < resolution:
				m.addConstr(x[lseg + lc + di] == 0.0)
				#print ("Set coverage of bp edge at index %d to 0." %(di))
		for srci in range(lsrc):
			cn_expr = gp.QuadExpr(0.0)
			cn_expr += w[0] * x[lseg + lc + ld + 2 * srci]
			cn_expr += w[0] * x[lseg + lc + ld + 2 * srci + 1]
			m.addQConstr(cn_expr <= remaining_CN[('src', srci)])
			
		# Occurrence of breakpoints in each cycle/path
		for seqi in range(lseg):
			m.addConstr(x[seqi] <= max_seq_repeat)

		# c: decomposition i is a cycle, and start at particular node
		c_names = []
		for ni in range(nnodes):
			c_names.append("c" + str(ni))
		c = m.addVars(nnodes, vtype = GRB.BINARY, name = c_names)

		# Relationship between c and x
		cycle_expr = gp.LinExpr(0.0)
		for ni in range(nnodes):
			cycle_expr += c[ni]
		for ni in range(len(endnode_list)):
			cycle_expr += x[lseg + lc + ld + 2 * lsrc + 2 * ni] # (s, v)
		for srci in range(lsrc): 
			cycle_expr += x[lseg + lc + ld + 2 * srci] # (s, v)
		m.addConstr(cycle_expr <= 1.0)

		# special request for c added for max_seq_repeat >= 2
		for node in g.nodes.keys():
			m.addConstr(c[node_order[node]] * x[g.nodes[node][0][0]] <= 1.0)
			
		# d: BFS/spanning tree order of the nodes in decomposition i
		d_names = []
		for ni in range(nnodes):
			d_names.append("d" + str(ni))
		d_names.append("ds")
		d_names.append("dt")
		d = m.addVars(nnodes + 2, lb = 0.0, ub = nnodes + 2, vtype = GRB.INTEGER, name = d_names) # including s, t, at the end
			
		# y: spanning tree indicator (directed)
		y1_names = []
		y2_names = []
		for ei in range(nedges):
			y1_names.append("y1-" + str(ei))
			y2_names.append("y2-" + str(ei))
		y1 = m.addVars(nedges, vtype = GRB.BINARY, name = y1_names) #small -> large
		y2 = m.addVars(nedges, vtype = GRB.BINARY, name = y2_names) #large -> small
			
		# Relationship between c and d
		for ni in range(nnodes):
			m.addConstr(d[ni] >= c[ni])
		d_expr = gp.LinExpr(0.0)
		for ni in range(nnodes):
			d_expr += c[ni]
		d_expr += d[nnodes]
		m.addConstr(d_expr <= 1.0)
			
		# Relationship between y and z:
		for ei in range(nedges):
			m.addConstr(y1[ei] <= z[0])
			m.addConstr(y2[ei] <= z[0])

		# Relationship between x, y and d
		for ei in range(nedges):
			m.addConstr(y1[ei] + y2[ei] <= x[ei])
		for di in range(ld):
			de = g.discordant_edges[di]
			if de[0] == de[3] and de[1] == de[4] and de[2] == de[5]: # exclude self loops
				m.addConstr(y1[lseg + lc + di] == 0)
				m.addConstr(y2[lseg + lc + di] == 0)
			
		t_expr_x = gp.LinExpr(0.0)
		t_expr_y = gp.LinExpr(0.0)
		t_expr_yd = gp.QuadExpr(0.0)
		for ni in range(len(endnode_list)):
			node = endnode_list[ni]
			t_expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * ni + 1]
			t_expr_y += y1[lseg + lc + ld + 2 * lsrc + 2 * ni + 1]
			t_expr_yd += y1[lseg + lc + ld + 2 * lsrc + 2 * ni + 1] * \
					(d[nnodes + 1] - d[node_order[node]]) # node -> t
			expr_x = gp.LinExpr(0.0)
			expr_y = gp.LinExpr(0.0)
			expr_xc = gp.QuadExpr(0.0)
			expr_yd = gp.QuadExpr(0.0)
			for seqi in g.nodes[node][0]:
				sseg = g.sequence_edges[seqi]
				node_ = (sseg[0], sseg[1], '-')
				if node_ == node:
					node_ = (sseg[0], sseg[2], '+')
				expr_x += x[seqi]
				expr_xc += x[seqi] * c[node_order[node]]
				if node_order[node_] <= node_order[node]:
					expr_y += y1[seqi]
					expr_yd += y1[seqi] * (d[node_order[node]] - d[node_order[node_]])
				else:
					expr_y += y2[seqi]
					expr_yd += y2[seqi] * (d[node_order[node]] - d[node_order[node_]])
						
			expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * ni] # from s
			expr_xc += x[lseg + lc + ld + 2 * lsrc + 2 * ni] * c[node_order[node]]
			expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * ni + 1] # to t
			expr_xc += x[lseg + lc + ld + 2 * lsrc + 2 * ni + 1] * c[node_order[node]]
			expr_y += y1[lseg + lc + ld + 2 * lsrc + 2 * ni] # from s
			expr_yd += y1[lseg + lc + ld + 2 * lsrc + 2 * ni] * \
								(d[node_order[node]] - d[nnodes])
			m.addConstr(expr_x * (nnodes + 2) >= d[node_order[node]])
			m.addConstr(expr_y <= 1.0)
			m.addConstr(expr_y * nedges + expr_xc >= expr_x)
			m.addConstr(expr_yd * nedges + expr_xc >= expr_x)
		for srci in range(lsrc):
			srce = g.source_edges[srci]
			t_expr_x += x[lseg + lc + ld + 2 * srci + 1]
			t_expr_y += y1[lseg + lc + ld + 2 * srci + 1]
			t_expr_yd += y1[lseg + lc + ld + 2 * srci + 1] * \
					(d[nnodes + 1] - d[node_order[(srce[3], srce[4], srce[5])]])
		m.addConstr(t_expr_x * (nnodes + 2) >= d[nnodes + 1])
		m.addConstr(t_expr_y <= 1.0)
		m.addConstr(t_expr_y * nedges >= t_expr_x)
		m.addConstr(t_expr_yd >= t_expr_x)
			
		for node in g.nodes.keys():
			if node not in endnode_list:
				expr_x = gp.LinExpr(0.0)
				expr_y = gp.LinExpr(0.0)
				expr_xc = gp.QuadExpr(0.0)
				expr_yd = gp.QuadExpr(0.0)
				for seqi in g.nodes[node][0]:
					sseg = g.sequence_edges[seqi]
					node_ = (sseg[0], sseg[1], '-')
					if node_ == node:
						node_ = (sseg[0], sseg[2], '+')
					expr_x += x[seqi]
					expr_xc += x[seqi] * c[node_order[node]]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[seqi]
						expr_yd += y1[seqi] * (d[node_order[node]] - d[node_order[node_]])
					else:
						expr_y += y2[seqi]
						expr_yd += y2[seqi] * (d[node_order[node]] - d[node_order[node_]])
				for ci in g.nodes[node][1]:
					ce = g.concordant_edges[ci]
					node_ = (ce[0], ce[1], ce[2])
					if node_ == node:
						node_ = (ce[3], ce[4], ce[5])
					expr_x += x[lseg + ci]
					expr_xc += x[lseg + ci] * c[node_order[node]]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[lseg + ci]
						expr_yd += y1[lseg + ci] * (d[node_order[node]] - d[node_order[node_]])
					else:
						expr_y += y2[lseg + ci]
						expr_yd += y2[lseg + ci] * (d[node_order[node]] - d[node_order[node_]])
				for di in g.nodes[node][2]:
					de = g.discordant_edges[di]
					node_ = (de[0], de[1], de[2])
					if node_ == node:
						node_ = (de[3], de[4], de[5])
					expr_x += x[lseg + lc + di]
					expr_xc += x[lseg + lc + di] * c[node_order[node]]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[lseg + lc + di]
						expr_yd += y1[lseg + lc + di] * (d[node_order[node]] - d[node_order[node_]])
					else:
						expr_y += y2[lseg + lc + di]
						expr_yd += y2[lseg + lc + di] * (d[node_order[node]] - d[node_order[node_]])
				for srci in g.nodes[node][3]:
					expr_x += x[lseg + lc + ld + 2 * srci]
					expr_x += x[lseg + lc + ld + 2 * srci + 1]
					expr_xc += x[lseg + lc + ld + 2 * srci] * c[node_order[node]]
					expr_xc += x[lseg + lc + ld + 2 * srci + 1] * c[node_order[node]]
					expr_y += y1[lseg + lc + ld + 2 * srci]
					expr_yd += y1[lseg + lc + ld + 2 * srci] * (d[node_order[node]] - d[nnodes])
				m.addConstr(expr_x * (nnodes + 2) >= d[node_order[node]])
				m.addConstr(expr_y <= 1.0)
				m.addConstr(expr_y * nedges + expr_xc >= expr_x)
				m.addConstr(expr_yd * nedges + expr_xc >= expr_x)

		# Subpath constraints
		for pi in range(len(pc_list)):
			path_constraint_ = pc_list[pi]
			for edge in path_constraint_.keys():
				if edge[0] == 's':
					m.addConstr(x[edge[1]] >= r[pi] * path_constraint_[edge])
				elif edge[0] == 'c':
					m.addConstr(x[lseg + edge[1]] >= r[pi] * path_constraint_[edge])
				else:
					m.addConstr(x[lseg + lc + edge[1]] >= r[pi] * path_constraint_[edge])
			
		m.setParam(GRB.Param.Threads, num_threads)
		m.setParam(GRB.Param.NonConvex, 2)
		m.setParam(GRB.Param.TimeLimit, 7200) # each breakpoint edge is assigned 5 minutes
		m.write(model_prefix + "_greedy_model" + str(cycle_id + 1) + "_alpha=" + str(alpha) + ".lp") 
		m.optimize()
		print('MS:', m.Status)
		if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
			print('Optimization was stopped with status %d' % m.Status)
			break
		else:
			cycle_id += 1
			total_weights_included = 0.0
			sol_z = m.getAttr('X', z)
			sol_w = m.getAttr('X', w)
			sol_d = m.getAttr('X', d)
			sol_r = m.getAttr('X', r)	
			if sol_z[0] >= 0.9:
				print ("Next cycle/Path exists; weight = %f" %(sol_w[0]))
				next_w = sol_w[0]
				if next_w < resolution:
					break
				sol_x = m.getAttr('X', x)
				sol_y1 = m.getAttr('X', y1)
				sol_y2 = m.getAttr('X', y2)
				sol_c = m.getAttr('X', c)
				cycle_flag = -1
				for ci in range(len(sol_c)):
					if sol_c[ci] >= 0.9:
						cycle_flag = ci
						break
				cycle = dict()
				for xi in range(len(sol_x)):
					cycle[('x', xi)] = sol_x[xi]
				for ci in range(len(sol_c)):
					cycle[('c', ci)] = sol_c[ci]
				for di in range(len(sol_d)):
					cycle[('d', di)] = sol_d[di]
				for yi in range(len(sol_y1)):
					cycle[('y1', yi)] = sol_y1[yi]
				for yi in range(len(sol_y2)):
					cycle[('y2', yi)] = sol_y2[yi]
				if cycle_flag == -1:
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if sol_x[xi] >= 0.9:
							x_xi = int(round(sol_x[xi]))
							if xi < lseg:
								print (cycle_id, 'path', 'seq', x_xi, g.sequence_edges[xi][:3])
								remaining_CN[('s', xi)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('s', xi)] < resolution:
									remaining_CN[('s', xi)] = 0.0
							elif xi < lseg + lc:
								print (cycle_id, 'path', 'concordant', x_xi, g.concordant_edges[xi - lseg][:6])
								remaining_CN[('c', xi - lseg)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('c', xi - lseg)] < resolution:
									remaining_CN[('c', xi - lseg)] = 0.0
							elif xi < lseg + lc + ld:
								print (cycle_id, 'path', 'discordant', x_xi, g.discordant_edges[xi - lseg - lc][:6])
								remaining_CN[('d', xi - lseg - lc)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('d', xi - lseg - lc)] < resolution:
									remaining_CN[('d', xi - lseg - lc)] = 0.0
							elif xi < lseg + lc + ld + 2 * lsrc:
								assert x_xi == 1
								if (xi - lseg - lc - ld) % 2 == 0:
									print (cycle_id, 'path', 'source', x_xi, g.source_edges[(xi - lseg - lc - ld) // 2][:6])
									remaining_CN[('src', (xi - lseg - lc - ld) // 2)] -= sol_w[0]
									if remaining_CN[('src', (xi - lseg - lc - ld) // 2)] < resolution:
										remaining_CN[('src', (xi - lseg - lc - ld) // 2)] = 0.0
								else:
									print (cycle_id, 'path', 'source', x_xi, g.source_edges[(xi - lseg - lc - ld - 1) // 2][:6])
									remaining_CN[('src', (xi - lseg - lc - ld - 1) // 2)] -= sol_w[0]
									if remaining_CN[('src', (xi - lseg - lc - ld - 1) // 2)] < resolution:
										remaining_CN[('src', (xi - lseg - lc - ld - 1) // 2)] = 0.0
							else:
								if (xi - lseg - lc - ld - 2 * lsrc) % 2 == 0:
									nsi = (xi - lseg - lc - ld - 2 * lsrc) // 2
									assert x_xi == 1
									print (cycle_id, 'path', 'source', x_xi, endnode_list[nsi])
								else:
									nti = (xi - lseg - lc - ld - 2 * lsrc - 1) // 2
									assert x_xi == 1
									print (cycle_id, 'path', 'source', x_xi, endnode_list[nti])
					for pi in range(len(pc_list)):
						if sol_r[pi] >= 0.9:
							path_constraints_s.append(pc_list[pi])
							unsatisfied_pc[pi] = -1
				else:
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if sol_x[xi] >= 0.9:
							x_xi = int(round(sol_x[xi]))
							if xi < lseg:
								print (cycle_id, 'cycle', 'seq', x_xi, g.sequence_edges[xi][:3])
								remaining_CN[('s', xi)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('s', xi)] < resolution:
									remaining_CN[('s', xi)] = 0.0
							elif xi < lseg + lc:
								print (cycle_id, 'cycle', 'concordant', x_xi, g.concordant_edges[xi - lseg][:6])
								remaining_CN[('c', xi - lseg)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('c', xi - lseg)] < resolution:
									remaining_CN[('c', xi - lseg)] = 0.0
							elif xi < lseg + lc + ld:
								print (cycle_id, 'cycle', 'discordant', x_xi, g.discordant_edges[xi - lseg - lc][:6])
								remaining_CN[('d', xi - lseg - lc)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('d', xi - lseg - lc)] < resolution:
									remaining_CN[('d', xi - lseg - lc)] = 0.0
							else:
								print ("Cyclic path cannot connect to source nodes.")
								os.abort()
					for pi in range(len(pc_list)):
						if sol_r[pi] >= 0.9:
							path_constraints_s.append(pc_list[pi])
							unsatisfied_pc[pi] = -1	
				for seqi in range(lseg):
					print (sol_x[seqi], sol_w[0], g.sequence_edges[seqi][-2])	
					total_weights_included += (sol_x[seqi] * sol_w[0] * g.sequence_edges[seqi][-2])
				print ("Total weights = ", total_weights_included, total_weights)
				remaining_weights -= total_weights_included
				if total_weights_included < 0.005 * total_weights:
					break
			else:
				break
		num_unsatisfied_pc = 0
		for i in range(len(pc_list)):
			if unsatisfied_pc[i] >= 0:
				num_unsatisfied_pc += 1
	#return remaining_weights


def maximize_weights_greedy_(g, total_weights, node_order, pc_list, alpha = 0.01,
			max_seq_repeat = 2, p_total_weight = 0.9, resolution = 0.1, num_threads = 16, model_prefix = ""):
	lseg = len(g.sequence_edges)
	lc = len(g.concordant_edges)
	ld = len(g.discordant_edges)
	lsrc = len(g.source_edges)
	print(lseg, lc, ld, lsrc)
	nedges = lseg + lc + ld + 2 * lsrc + 2 * len(g.endnodes)
	nnodes = len(g.nodes)
	endnode_list = [node for node in g.endnodes.keys()]

	remaining_weights = total_weights
	unsatisfied_pc = [i for i in range(len(pc_list))]
	remaining_CN = dict()
	for seqi in range(lseg):
		remaining_CN[('s', seqi)] = g.sequence_edges[seqi][-1]
	for ci in range(lc):
		remaining_CN[('c', ci)] = g.concordant_edges[ci][-1]
	for di in range(ld):
		remaining_CN[('d', di)] = g.discordant_edges[di][-1]
	for srci in range(lsrc):
		remaining_CN[('src', srci)] = g.source_edges[srci][-1]

	next_w = resolution * 1.1
	cycle_id = 0
	num_unsatisfied_pc = len(pc_list)
	while next_w >= resolution and (remaining_weights > (1.0 - p_total_weight) * total_weights or num_unsatisfied_pc > math.floor(0.1 * len(pc_list))):
		pp = 1.0
		if alpha > 0 and num_unsatisfied_pc > 0:
			pp = alpha * remaining_weights / num_unsatisfied_pc # multi - objective optimization parameter
		print ("Cycle id = ", cycle_id)
		print ("Remaining weights = ", remaining_weights, total_weights)
		print ("Num unsatisfied path constraints = ", num_unsatisfied_pc, len(pc_list))
		print ("Path constraints factor = ", pp)

		# Gurobi model
		m = gp.Model(model_prefix + "_cycle_decomposition_greedy_" + str(cycle_id + 1))

		# z[i]: indicating whether cycle or path i exists
		z = m.addVars(1, vtype = GRB.BINARY, name = ["z0"])

		# w[i]: the weight of cycle or path i, continuous variable
		w = m.addVars(1, lb = 0.0, ub = g.max_cn, vtype = GRB.CONTINUOUS, name = ["w0"])

		# Relationship between w[i] and z[i]
		m.addConstr(w[0] <= z[0] * g.max_cn)
		m.addConstr(w[0] >= z[0] * resolution)

		# x: the number of times an edge occur in cycle or path i
		x_names = []
		for ei in range(nedges):
			x_names.append("x" + str(ei))
		x = m.addVars(nedges, lb = 0.0, ub = 10.0, vtype = GRB.INTEGER, name = x_names)

		# Subpath constraints
		r_names = []
		for pi in range(len(pc_list)):
			r_names.append("r" + str(pi))
		r = m.addVars(len(pc_list), vtype = GRB.BINARY, name = r_names)
			
		# Objective: maximize the total weight + num subpath constraints satisfied
		obj = gp.QuadExpr(0.0)
		obj += w[0]
		m.setObjective(obj, GRB.MAXIMIZE)

		# Eulerian constraint
		for node in g.nodes.keys():
			if node in endnode_list:
				m.addConstr(x[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)] + \
						x[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node) + 1] \
						== x[g.nodes[node][0][0]])
			else:
				ec_expr = gp.LinExpr(0.0)
				for seqi in g.nodes[node][0]:
					ec_expr += x[seqi]
				for ci in g.nodes[node][1]:
					ec_expr -= x[lseg + ci]
				for di in g.nodes[node][2]:
					de = g.discordant_edges[di]
					if de[0] == de[3] and de[1] == de[4] and de[2] == de[5]:
						ec_expr -= 2 * x[lseg + lc + di]
					else:
						ec_expr -= x[lseg + lc + di]
				for srci in g.nodes[node][3]:
					ec_expr -= x[lseg + lc + ld + 2 * srci] # connected to s
					ec_expr -= x[lseg + lc + ld + 2 * srci + 1] # connected to t
				m.addConstr(ec_expr == 0.0)
		path_expr = gp.LinExpr(0.0)
		for ni in range(len(endnode_list)):
			path_expr += x[lseg + lc + ld + 2 * lsrc + 2 * ni] # (s, v)
			path_expr -= x[lseg + lc + ld + 2 * lsrc + 2 * ni + 1] # (v, t)
		for srci in range(lsrc): 
			path_expr += x[lseg + lc + ld + 2 * srci] # (s, v)
			path_expr -= x[lseg + lc + ld + 2 * srci + 1] # (v, t)
		m.addConstr(path_expr == 0.0)

		# CN constraint
		length_expr = gp.LinExpr(0.0)
		for seqi in range(lseg):
			m.addQConstr(w[0] * x[seqi] <= remaining_CN[('s', seqi)])
			length_expr += x[seqi] * g.sequence_edges[seqi][-2]
		m.addConstr(length_expr >= 10000.0)
		for ci in range(lc):
			m.addQConstr(w[0] * x[lseg + ci] <= remaining_CN[('c', ci)])
		for di in range(ld):
			m.addQConstr(w[0] * x[lseg + lc + di] <= remaining_CN[('d', di)])
			if remaining_CN[('d', di)] < resolution:
				m.addConstr(x[lseg + lc + di] == 0.0)
				#print ("Set coverage of bp edge at index %d to 0." %(di))
		for srci in range(lsrc):
			cn_expr = gp.QuadExpr(0.0)
			cn_expr += w[0] * x[lseg + lc + ld + 2 * srci]
			cn_expr += w[0] * x[lseg + lc + ld + 2 * srci + 1]
			m.addQConstr(cn_expr <= remaining_CN[('src', srci)])
			
		# Occurrence of breakpoints in each cycle/path
		for seqi in range(lseg):
			m.addConstr(x[seqi] <= max_seq_repeat)

		# c: decomposition i is a cycle, and start at particular node
		c_names = []
		for ni in range(nnodes):
			c_names.append("c" + str(ni))
		c = m.addVars(nnodes, vtype = GRB.BINARY, name = c_names)

		# Relationship between c and x
		cycle_expr = gp.LinExpr(0.0)
		for ni in range(nnodes):
			cycle_expr += c[ni]
		for ni in range(len(endnode_list)):
			cycle_expr += x[lseg + lc + ld + 2 * lsrc + 2 * ni] # (s, v)
		for srci in range(lsrc): 
			cycle_expr += x[lseg + lc + ld + 2 * srci] # (s, v)
		m.addConstr(cycle_expr <= 1.0)

		# special request for c added for max_seq_repeat >= 2
		for node in g.nodes.keys():
			m.addConstr(c[node_order[node]] * x[g.nodes[node][0][0]] <= 1.0)
			
		# d: BFS/spanning tree order of the nodes in decomposition i
		d_names = []
		for ni in range(nnodes):
			d_names.append("d" + str(ni))
		d_names.append("ds")
		d_names.append("dt")
		d = m.addVars(nnodes + 2, lb = 0.0, ub = nnodes + 2, vtype = GRB.INTEGER, name = d_names) # including s, t, at the end
			
		# y: spanning tree indicator (directed)
		y1_names = []
		y2_names = []
		for ei in range(nedges):
			y1_names.append("y1-" + str(ei))
			y2_names.append("y2-" + str(ei))
		y1 = m.addVars(nedges, vtype = GRB.BINARY, name = y1_names) #small -> large
		y2 = m.addVars(nedges, vtype = GRB.BINARY, name = y2_names) #large -> small
			
		# Relationship between c and d
		for ni in range(nnodes):
			m.addConstr(d[ni] >= c[ni])
		d_expr = gp.LinExpr(0.0)
		for ni in range(nnodes):
			d_expr += c[ni]
		d_expr += d[nnodes]
		m.addConstr(d_expr <= 1.0)
			
		# Relationship between y and z:
		for ei in range(nedges):
			m.addConstr(y1[ei] <= z[0])
			m.addConstr(y2[ei] <= z[0])

		# Relationship between x, y and d
		for ei in range(nedges):
			m.addConstr(y1[ei] + y2[ei] <= x[ei])
		for di in range(ld):
			de = g.discordant_edges[di]
			if de[0] == de[3] and de[1] == de[4] and de[2] == de[5]: # exclude self loops
				m.addConstr(y1[lseg + lc + di] == 0)
				m.addConstr(y2[lseg + lc + di] == 0)
			
		t_expr_x = gp.LinExpr(0.0)
		t_expr_y = gp.LinExpr(0.0)
		t_expr_yd = gp.QuadExpr(0.0)
		for ni in range(len(endnode_list)):
			node = endnode_list[ni]
			t_expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * ni + 1]
			t_expr_y += y1[lseg + lc + ld + 2 * lsrc + 2 * ni + 1]
			t_expr_yd += y1[lseg + lc + ld + 2 * lsrc + 2 * ni + 1] * \
					(d[nnodes + 1] - d[node_order[node]]) # node -> t
			expr_x = gp.LinExpr(0.0)
			expr_y = gp.LinExpr(0.0)
			expr_xc = gp.QuadExpr(0.0)
			expr_yd = gp.QuadExpr(0.0)
			for seqi in g.nodes[node][0]:
				sseg = g.sequence_edges[seqi]
				node_ = (sseg[0], sseg[1], '-')
				if node_ == node:
					node_ = (sseg[0], sseg[2], '+')
				expr_x += x[seqi]
				expr_xc += x[seqi] * c[node_order[node]]
				if node_order[node_] <= node_order[node]:
					expr_y += y1[seqi]
					expr_yd += y1[seqi] * (d[node_order[node]] - d[node_order[node_]])
				else:
					expr_y += y2[seqi]
					expr_yd += y2[seqi] * (d[node_order[node]] - d[node_order[node_]])
						
			expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * ni] # from s
			expr_xc += x[lseg + lc + ld + 2 * lsrc + 2 * ni] * c[node_order[node]]
			expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * ni + 1] # to t
			expr_xc += x[lseg + lc + ld + 2 * lsrc + 2 * ni + 1] * c[node_order[node]]
			expr_y += y1[lseg + lc + ld + 2 * lsrc + 2 * ni] # from s
			expr_yd += y1[lseg + lc + ld + 2 * lsrc + 2 * ni] * \
								(d[node_order[node]] - d[nnodes])
			m.addConstr(expr_x * (nnodes + 2) >= d[node_order[node]])
			m.addConstr(expr_y <= 1.0)
			m.addConstr(expr_y * nedges + expr_xc >= expr_x)
			m.addConstr(expr_yd * nedges + expr_xc >= expr_x)
		for srci in range(lsrc):
			srce = g.source_edges[srci]
			t_expr_x += x[lseg + lc + ld + 2 * srci + 1]
			t_expr_y += y1[lseg + lc + ld + 2 * srci + 1]
			t_expr_yd += y1[lseg + lc + ld + 2 * srci + 1] * \
					(d[nnodes + 1] - d[node_order[(srce[3], srce[4], srce[5])]])
		m.addConstr(t_expr_x * (nnodes + 2) >= d[nnodes + 1])
		m.addConstr(t_expr_y <= 1.0)
		m.addConstr(t_expr_y * nedges >= t_expr_x)
		m.addConstr(t_expr_yd >= t_expr_x)
			
		for node in g.nodes.keys():
			if node not in endnode_list:
				expr_x = gp.LinExpr(0.0)
				expr_y = gp.LinExpr(0.0)
				expr_xc = gp.QuadExpr(0.0)
				expr_yd = gp.QuadExpr(0.0)
				for seqi in g.nodes[node][0]:
					sseg = g.sequence_edges[seqi]
					node_ = (sseg[0], sseg[1], '-')
					if node_ == node:
						node_ = (sseg[0], sseg[2], '+')
					expr_x += x[seqi]
					expr_xc += x[seqi] * c[node_order[node]]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[seqi]
						expr_yd += y1[seqi] * (d[node_order[node]] - d[node_order[node_]])
					else:
						expr_y += y2[seqi]
						expr_yd += y2[seqi] * (d[node_order[node]] - d[node_order[node_]])
				for ci in g.nodes[node][1]:
					ce = g.concordant_edges[ci]
					node_ = (ce[0], ce[1], ce[2])
					if node_ == node:
						node_ = (ce[3], ce[4], ce[5])
					expr_x += x[lseg + ci]
					expr_xc += x[lseg + ci] * c[node_order[node]]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[lseg + ci]
						expr_yd += y1[lseg + ci] * (d[node_order[node]] - d[node_order[node_]])
					else:
						expr_y += y2[lseg + ci]
						expr_yd += y2[lseg + ci] * (d[node_order[node]] - d[node_order[node_]])
				for di in g.nodes[node][2]:
					de = g.discordant_edges[di]
					node_ = (de[0], de[1], de[2])
					if node_ == node:
						node_ = (de[3], de[4], de[5])
					expr_x += x[lseg + lc + di]
					expr_xc += x[lseg + lc + di] * c[node_order[node]]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[lseg + lc + di]
						expr_yd += y1[lseg + lc + di] * (d[node_order[node]] - d[node_order[node_]])
					else:
						expr_y += y2[lseg + lc + di]
						expr_yd += y2[lseg + lc + di] * (d[node_order[node]] - d[node_order[node_]])
				for srci in g.nodes[node][3]:
					expr_x += x[lseg + lc + ld + 2 * srci]
					expr_x += x[lseg + lc + ld + 2 * srci + 1]
					expr_xc += x[lseg + lc + ld + 2 * srci] * c[node_order[node]]
					expr_xc += x[lseg + lc + ld + 2 * srci + 1] * c[node_order[node]]
					expr_y += y1[lseg + lc + ld + 2 * srci]
					expr_yd += y1[lseg + lc + ld + 2 * srci] * (d[node_order[node]] - d[nnodes])
				m.addConstr(expr_x * (nnodes + 2) >= d[node_order[node]])
				m.addConstr(expr_y <= 1.0)
				m.addConstr(expr_y * nedges + expr_xc >= expr_x)
				m.addConstr(expr_yd * nedges + expr_xc >= expr_x)

		# Subpath constraints
		for pi in range(len(pc_list)):
			path_constraint_ = pc_list[pi]
			for edge in path_constraint_.keys():
				if edge[0] == 's':
					m.addConstr(x[edge[1]] >= r[pi] * path_constraint_[edge])
				elif edge[0] == 'c':
					m.addConstr(x[lseg + edge[1]] >= r[pi] * path_constraint_[edge])
				else:
					m.addConstr(x[lseg + lc + edge[1]] >= r[pi] * path_constraint_[edge])
			
		m.setParam(GRB.Param.Threads, num_threads)
		m.setParam(GRB.Param.NonConvex, 2)
		m.setParam(GRB.Param.TimeLimit, 7200 / num_threads) # each breakpoint edge is assigned 5 minutes
		m.write(model_prefix + "_greedy_model" + str(cycle_id + 1) + "_alpha=" + str(alpha) + ".lp") 
		m.optimize()
		print('MS:', m.Status)
		if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
			print('Optimization was stopped with status %d' % m.Status)
			break
		else:
			cycle_id += 1
			total_weights_included = 0.0
			sol_z = m.getAttr('X', z)
			sol_w = m.getAttr('X', w)
			sol_d = m.getAttr('X', d)
			sol_r = m.getAttr('X', r)	
			if sol_z[0] >= 0.9:
				print ("Next cycle/Path exists; weight = %f" %(sol_w[0]))
				next_w = sol_w[0]
				if next_w < resolution:
					break
				sol_x = m.getAttr('X', x)
				sol_y1 = m.getAttr('X', y1)
				sol_y2 = m.getAttr('X', y2)
				sol_c = m.getAttr('X', c)
				cycle_flag = -1
				for ci in range(len(sol_c)):
					if sol_c[ci] >= 0.9:
						cycle_flag = ci
						break
				cycle = dict()
				for xi in range(len(sol_x)):
					cycle[('x', xi)] = sol_x[xi]
				for ci in range(len(sol_c)):
					cycle[('c', ci)] = sol_c[ci]
				for di in range(len(sol_d)):
					cycle[('d', di)] = sol_d[di]
				for yi in range(len(sol_y1)):
					cycle[('y1', yi)] = sol_y1[yi]
				for yi in range(len(sol_y2)):
					cycle[('y2', yi)] = sol_y2[yi]
				if cycle_flag == -1:
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if sol_x[xi] >= 0.9:
							x_xi = int(round(sol_x[xi]))
							if xi < lseg:
								print (cycle_id, 'path', 'seq', x_xi, g.sequence_edges[xi][:3])
								remaining_CN[('s', xi)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('s', xi)] < resolution:
									remaining_CN[('s', xi)] = 0.0
							elif xi < lseg + lc:
								print (cycle_id, 'path', 'concordant', x_xi, g.concordant_edges[xi - lseg][:6])
								remaining_CN[('c', xi - lseg)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('c', xi - lseg)] < resolution:
									remaining_CN[('c', xi - lseg)] = 0.0
							elif xi < lseg + lc + ld:
								print (cycle_id, 'path', 'discordant', x_xi, g.discordant_edges[xi - lseg - lc][:6])
								remaining_CN[('d', xi - lseg - lc)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('d', xi - lseg - lc)] < resolution:
									remaining_CN[('d', xi - lseg - lc)] = 0.0
							elif xi < lseg + lc + ld + 2 * lsrc:
								assert x_xi == 1
								if (xi - lseg - lc - ld) % 2 == 0:
									print (cycle_id, 'path', 'source', x_xi, g.source_edges[(xi - lseg - lc - ld) // 2][:6])
									remaining_CN[('src', (xi - lseg - lc - ld) // 2)] -= sol_w[0]
									if remaining_CN[('src', (xi - lseg - lc - ld) // 2)] < resolution:
										remaining_CN[('src', (xi - lseg - lc - ld) // 2)] = 0.0
								else:
									print (cycle_id, 'path', 'source', x_xi, g.source_edges[(xi - lseg - lc - ld - 1) // 2][:6])
									remaining_CN[('src', (xi - lseg - lc - ld - 1) // 2)] -= sol_w[0]
									if remaining_CN[('src', (xi - lseg - lc - ld - 1) // 2)] < resolution:
										remaining_CN[('src', (xi - lseg - lc - ld - 1) // 2)] = 0.0
							else:
								if (xi - lseg - lc - ld - 2 * lsrc) % 2 == 0:
									nsi = (xi - lseg - lc - ld - 2 * lsrc) // 2
									assert x_xi == 1
									print (cycle_id, 'path', 'source', x_xi, endnode_list[nsi])
								else:
									nti = (xi - lseg - lc - ld - 2 * lsrc - 1) // 2
									assert x_xi == 1
									print (cycle_id, 'path', 'source', x_xi, endnode_list[nti])
					for pi in range(len(pc_list)):
						if sol_r[pi] >= 0.9:
							path_constraints_s.append(pc_list[pi])
							unsatisfied_pc[pi] = -1
				else:
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if sol_x[xi] >= 0.9:
							x_xi = int(round(sol_x[xi]))
							if xi < lseg:
								print (cycle_id, 'cycle', 'seq', x_xi, g.sequence_edges[xi][:3])
								remaining_CN[('s', xi)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('s', xi)] < resolution:
									remaining_CN[('s', xi)] = 0.0
							elif xi < lseg + lc:
								print (cycle_id, 'cycle', 'concordant', x_xi, g.concordant_edges[xi - lseg][:6])
								remaining_CN[('c', xi - lseg)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('c', xi - lseg)] < resolution:
									remaining_CN[('c', xi - lseg)] = 0.0
							elif xi < lseg + lc + ld:
								print (cycle_id, 'cycle', 'discordant', x_xi, g.discordant_edges[xi - lseg - lc][:6])
								remaining_CN[('d', xi - lseg - lc)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('d', xi - lseg - lc)] < resolution:
									remaining_CN[('d', xi - lseg - lc)] = 0.0
							else:
								print ("Cyclic path cannot connect to source nodes.")
								os.abort()
					for pi in range(len(pc_list)):
						if sol_r[pi] >= 0.9:
							path_constraints_s.append(pc_list[pi])
							unsatisfied_pc[pi] = -1	
				for seqi in range(lseg):
					print (sol_x[seqi], sol_w[0], g.sequence_edges[seqi][-2])	
					total_weights_included += (sol_x[seqi] * sol_w[0] * g.sequence_edges[seqi][-2])
				print ("Total weights = ", total_weights_included, total_weights)
				remaining_weights -= total_weights_included
				if total_weights_included < 0.005 * total_weights:
					break
			else:
				break
		num_unsatisfied_pc = 0
		for i in range(len(pc_list)):
			if unsatisfied_pc[i] >= 0:
				num_unsatisfied_pc += 1
	#return remaining_weights

