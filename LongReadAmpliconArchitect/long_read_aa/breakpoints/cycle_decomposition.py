"""
Functions used for cycle decomposition
"""
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
		node_ = (sseg[0], sseg[1], '-')
		if node_ == node:
			node_ = (sseg[0], sseg[2], '+')
		yi = endnode_list.index(node)
		if node_order[node] < node_order[node_]:
			m.addConstr(y2[seqi] + y1[lseg + lc + ld + lsrc + yi] == y1[seqi] + y2[lseg + lc + ld + lsrc + yi])
		else:
			m.addConstr(y1[seqi] + y1[lseg + lc + ld + lsrc + yi] == y2[seqi] + y2[lseg + lc + ld + lsrc + yi])
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
	residual_network = copy.deepcopy(g)
	if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
		print('Optimization was stopped with status %d' % m.Status)
		return GRB.INFEASIBLE, flow_network, residual_network
	else:
		sol_y1 = m.getAttr('X', y1)
		sol_y2 = m.getAttr('X', y2)
		for ni in range(len(endnode_list)):
			node = endnode_list[ni]
			sseg = g.sequence_edges[g.nodes[node][0][0]]
			if sol_y1[lseg + lc + ld + lsrc + ni] + sol_y2[lseg + lc + ld + lsrc + ni] > 0:
				flow_network.add_endnode(node)
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
					print (de, sol_y1[yi], sol_y2[yi])
					flow_network.add_discordant_edge(de[0], de[1], de[2], de[3], de[4], de[5], cn = (sol_y1[yi] + sol_y2[yi]))
				elif yi < lseg + lc + ld + lsrc:
					srce = g.source_edges[yi - lseg - lc - ld]
					residual_network.source_edges[yi - lseg - lc - ld][-1] -= (sol_y1[yi] + sol_y2[yi])
					flow_network.add_source_edge(srce[3], srce[4], srce[5], cn = (sol_y1[yi] + sol_y2[yi]))
		return m.Status, flow_network, residual_network		


#def minimize_cycles_max_flow(g, )
	"""
	Given a breakpoint graph, compute the maximum s-t flow
	
	"""
#	return status, resedual, cycles


