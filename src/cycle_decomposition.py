"""
Functions used for cycle decomposition
"""
import math
import random
import gurobipy as gp
from gurobipy import GRB

from breakpoint_graph import *
import global_names


def minimize_cycles(amplicon_id, g, k, total_weights, node_order, pc_list, max_seq_repeat = 2,
		p_total_weight = 0.9, p_bp_cn = 0.9, num_threads = -1, time_limit = 7200, model_prefix = ""):
	"""
	Cycle decomposition by minimizing the number of cycles/paths

	amplicon_id: integer, amplicon ID
	g: breakpoint graph (object)
	k: integer, maximum mumber of cycles/paths allowed in cycle decomposition
	total_weights: float, total length-weighted CN in breakpoint graph g
	node_order: dict maps each node in the input breakpoint graphg to a distinct integer, indicating a total order of the nodes in g
	pc_list: list of subpath constraints to be satisfied, each as a dict that maps an edge to its multiplicity
		*** note that all subpath constraints in this list are required to be satisfied ***
		*** otherwise will return infeasible ***
	max_seq_repeat: integer, maximum multiplicity (num occurrence) allowed for each sequence edge in each cycle/path, 
			default value is 2
	p_total_weight: float between (0, 1), minimum proportion of length-weighted CN to be covered by the resulting cycles or paths, 
			default value is 0.9
	p_bp_cn: float float between (0, 1), minimum proportion of CN for each discordant edge to be covered by the resulting cycles or paths, 
			default value is 0.9
	num_threads: integer, number of working threads for gurobipy, by default it tries to use up all available cores 
	time_limit: integer, maximum allowed running time, in seconds, default is 7200 (2 hour)
	model_prefix: output prefix for gurobi *.lp model

	Returns: (1) Status of gurobi optimization model (usually 2 - optimal; 3 - infeasible; 9 - suboptimal/reached time limit)
		(2) Total length weighted CN in resulting cycles/paths
		(3) Total num subpath constraints satisfied by resulting cycles/paths
		(4) List of cycles, each as a dict which maps an edge to its multiplicity in the cycle
		(5) List of the corresponding CN of the above cycles
		(6) Subpath constraints (indices) satisfied by each cycle
		(7) List of paths, each as a dict which maps an edge to its multiplicity in the path
		(8) List of the corresponding CN of the above paths
		(9) Subpath constraints (indices) satisfied by each path
	"""
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Regular cycle decomposition with at most %d cycles/paths allowed." %k)
	lseg = len(g.sequence_edges)
	lc = len(g.concordant_edges)
	ld = len(g.discordant_edges)
	lsrc = len(g.source_edges)
	nnodes = len(g.nodes)
	nedges = lseg + lc + ld + 2 * lsrc + 2 * len(g.endnodes)
	endnode_list = [node for node in g.endnodes.keys()]
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tNum nodes to be used in QP = %d." %nnodes)
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tNum edges to be used in QP = %d." %nedges)

	# Gurobi model
	m = gp.Model(model_prefix + "_amplicon" + str(amplicon_id) + "_cycle_decomposition_k=" + str(k))
		
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
	obj = gp.QuadExpr(1.0)
	for i in range(k):
		obj += z[i]
		for seqi in range(lseg):
			obj -= (x[seqi * k + i] * w[i] * g.sequence_edges[seqi][-2] / total_weights)
	m.setObjective(obj, GRB.MINIMIZE)

	# Must include at least 0.9 * total CN weights
	total_weights_expr = gp.QuadExpr(0.0)
	for i in range(k):
		for seqi in range(lseg):
			total_weights_expr += (x[seqi * k + i] * w[i] * g.sequence_edges[seqi][-2])
	m.addConstr(total_weights_expr >= p_total_weight * total_weights)

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
					ec_expr -= x[(lseg + lc + di) * k + i]
				for srci in g.nodes[node][3]:
					ec_expr -= x[(lseg + lc + ld + 2 * srci) * k + i] # connected to s
					ec_expr -= x[(lseg + lc + ld + 2 * srci) * k + k + i] # connected to t
				m.addConstr(ec_expr == 0.0)
	for i in range(k):
		path_expr = gp.LinExpr(0.0)
		for enodei in range(len(endnode_list)):
			path_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * enodei) * k + i] # (s, v)
			path_expr -= x[(lseg + lc + ld + 2 * lsrc + 2 * enodei) * k + k + i] # (v, t)
		for srci in range(lsrc): 
			path_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
			path_expr -= x[(lseg + lc + ld + 2 * srci) * k + k + i] # (v, t)
		m.addConstr(path_expr == 0.0)

	# CN constraint
	for seqi in range(lseg):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[seqi * k + i]
		m.addQConstr(cn_expr <= g.sequence_edges[seqi][-1])
	for ci in range(lc):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[(lseg + ci) * k + i]
		m.addQConstr(cn_expr <= g.concordant_edges[ci][-1])
	for di in range(ld):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[(lseg + lc + di) * k + i]
		m.addQConstr(cn_expr <= g.discordant_edges[di][-1])
		m.addQConstr(cn_expr >= p_bp_cn * g.discordant_edges[di][-1])
	for srci in range(lsrc):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + i]
			cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + k + i]
		m.addQConstr(cn_expr <= g.source_edges[srci][-1])
			
	# Occurrence of breakpoints in each cycle/path
	for i in range(k):
		for seqi in range(lseg):
			m.addConstr(x[seqi * k + i] <= max_seq_repeat)
	for i in range(k):
		sum_inv = gp.LinExpr(1.0)
		for di in range(ld):
			dedge = g.discordant_edges[di]
			if dedge[2] == dedge[5]:
				sum_inv += x[(lseg + lc + di) * k + i]
		for di in range(ld):
			m.addConstr(x[(lseg + lc + di) * k + i] <= 2)
			m.addConstr(x[(lseg + lc + di) * k + i] <= sum_inv)

	# c: decomposition i is a cycle, and start at particular node
	c_names = []
	for ni in range(nnodes):
		for i in range(k):
			c_names.append("c" + str(ni) + "," + str(i))
	c = m.addVars(k * nnodes, vtype = GRB.BINARY, name = c_names)

	# Relationship between c and x
	for i in range(k):
		cycle_expr = gp.LinExpr(0.0)
		for ni in range(nnodes):
			cycle_expr += c[ni * k + i]
		for eni in range(len(endnode_list)):
			cycle_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * eni) * k + i] # (s, v)
		for srci in range(lsrc): 
			cycle_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
		m.addConstr(cycle_expr <= 1.0)

	# special request for c added for max_seq_repeat >= 2
	if max_seq_repeat >= 2:
		for i in range(k):
			expr_xc = gp.QuadExpr(0.0)
			for node in g.nodes.keys():
				for ci in set(g.nodes[node][1]):
					expr_xc += c[k * node_order[node] + i] * x[(lseg + ci) * k + i]
				for di in set(g.nodes[node][2]):	
					expr_xc += c[k * node_order[node] + i] * x[(lseg + lc + di) * k + i]
			m.addConstr(expr_xc <= 1.0)
			
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
			
	# Relationship between c and d
	for i in range(k):
		for ni in range(nnodes):
			m.addConstr(d[ni * k + i] >= c[ni * k + i])
	for i in range(k):
		d_expr = gp.LinExpr(0.0)
		for ni in range(nnodes):
			d_expr += c[ni * k + i]
		d_expr += d[k * nnodes + i] # d_s,i
		m.addConstr(d_expr <= 1.0)
			
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
			dedge = g.discordant_edges[di]
			if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
				m.addConstr(y1[(lseg + lc + di) * k + i] == 0)
				m.addConstr(y2[(lseg + lc + di) * k + i] == 0)
	for i in range(k):
		t_expr_x = gp.LinExpr(0.0)
		t_expr_y = gp.LinExpr(0.0)
		t_expr_yd = gp.QuadExpr(0.0)
		for node in endnode_list:
			t_expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i]
			t_expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i]
			t_expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i] * \
					(d[k * nnodes + k + i] - d[k * node_order[node] + i]) # node -> t
			expr_x = gp.LinExpr(0.0)
			expr_y = gp.LinExpr(0.0)
			expr_xc = gp.QuadExpr(0.0)
			expr_yd = gp.QuadExpr(0.0)
			for seqi in g.nodes[node][0]:
				sseg = g.sequence_edges[seqi]
				node_ = (sseg[0], sseg[1], '-')
				if node_ == node:
					node_ = (sseg[0], sseg[2], '+')
				expr_x += x[seqi * k + i]
				expr_xc += x[seqi * k + i] * c[k * node_order[node] + i]
				if node_order[node_] <= node_order[node]:
					expr_y += y1[seqi * k + i]
					expr_yd += y1[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				else:
					expr_y += y2[seqi * k + i]
					expr_yd += y2[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						
			expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + i] # from s
			expr_xc += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + i] * c[k * node_order[node] + i]
			expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i] # to t
			expr_xc += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i] * c[k * node_order[node] + i]
			expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + i] # from s
			expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + i] * \
					(d[k * node_order[node] + i] - d[k * nnodes + i])
			m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
			m.addConstr(expr_y <= 1.0)
			m.addConstr(expr_y * nedges * k + expr_xc >= expr_x)
			m.addConstr(expr_yd * nedges * k + expr_xc >= expr_x)
			
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
				expr_xc = gp.QuadExpr(0.0)
				expr_yd = gp.QuadExpr(0.0)
				for seqi in g.nodes[node][0]:
					sseg = g.sequence_edges[seqi]
					node_ = (sseg[0], sseg[1], '-')
					if node_ == node:
						node_ = (sseg[0], sseg[2], '+')
					expr_x += x[seqi * k + i]
					expr_xc += x[seqi * k + i] * c[k * node_order[node] + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[seqi * k + i]
						expr_yd += y1[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[seqi * k + i]
						expr_yd += y2[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				for ci in g.nodes[node][1]:
					cedge = g.concordant_edges[ci]
					node_ = (cedge[0], cedge[1], cedge[2])
					if node_ == node:
						node_ = (cedge[3], cedge[4], cedge[5])
					expr_x += x[(lseg + ci) * k + i]
					expr_xc += x[(lseg + ci) * k + i] * c[k * node_order[node] + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[(lseg + ci) * k + i]
						expr_yd += y1[(lseg + ci) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[(lseg + ci) * k + i]
						expr_yd += y2[(lseg + ci) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				for di in g.nodes[node][2]:
					dedge = g.discordant_edges[di]
					node_ = (dedge[0], dedge[1], dedge[2])
					if node_ == node:
						node_ = (dedge[3], dedge[4], dedge[5])
					expr_x += x[(lseg + lc + di) * k + i]
					expr_xc += x[(lseg + lc + di) * k + i] * c[k * node_order[node] + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[(lseg + lc + di) * k + i]
						expr_yd += y1[(lseg + lc + di) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[(lseg + lc + di) * k + i]
						expr_yd += y2[(lseg + lc + di) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				for srci in g.nodes[node][3]:
					expr_x += x[(lseg + lc + ld + 2 * srci) * k + i]
					expr_x += x[(lseg + lc + ld + 2 * srci) * k + k + i]
					expr_xc += x[(lseg + lc + ld + 2 * srci) * k + i] * c[k * node_order[node] + i]
					expr_xc += x[(lseg + lc + ld + 2 * srci) * k + k + i] * c[k * node_order[node] + i]
					expr_y += y1[(lseg + lc + ld + 2 * srci) * k + i]
					expr_yd += y1[(lseg + lc + ld + 2 * srci) * k + i] * \
							(d[k * node_order[node] + i] - d[k * nnodes + i])
				m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
				m.addConstr(expr_y <= 1.0)
				m.addConstr(expr_y * nedges * k + expr_xc >= expr_x)
				m.addConstr(expr_yd * nedges * k + expr_xc >= expr_x)

	# Subpath constraints
	r_names = []
	for pi in range(len(pc_list)):
		for i in range(k):
			r_names.append("r" + str(pi) + "," + str(i))
	r = m.addVars(k * len(pc_list), vtype = GRB.BINARY, name = r_names)
	for pi in range(len(pc_list)):
		path_constraint_ = pc_list[pi]
		sum_r = gp.LinExpr(0.0)
		for ri in range(pi * k, (pi + 1) * k):
			sum_r += r[ri]
		m.addConstr(sum_r >= 1.0)	
		for edge in path_constraint_.keys():
			for i in range(k):
				if edge[0] == 's':
					m.addConstr(x[edge[1] * k + i] >= r[pi * k + i] * path_constraint_[edge])
				elif edge[0] == 'c':
					m.addConstr(x[(lseg + edge[1]) * k + i] >= r[pi * k + i] * path_constraint_[edge])
				else:
					m.addConstr(x[(lseg + lc + edge[1]) * k + i] >= r[pi * k + i] * path_constraint_[edge])

	m.setParam(GRB.Param.LogToConsole, 0)
	if num_threads > 0:
		m.setParam(GRB.Param.Threads, num_threads)
	m.setParam(GRB.Param.NonConvex, 2)
	m.setParam(GRB.Param.TimeLimit, max(time_limit, ld * 300)) # each breakpoint edge is assigned 5 minutes
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCompleted gurobi model setup.")
	lp_fn = model_prefix + "_amplicon" + str(amplicon_id) + "_model.lp"
	m.write(lp_fn)
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tWrote model to file: %s." %lp_fn)
	log_fn = lp_fn[:-2] + "log"
	m.setParam(GRB.Param.LogFile, log_fn)	
	m.optimize()
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCompleted optimization with status %d." %m.Status)
	
	if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tModel is infeasible.")
		return GRB.INFEASIBLE, 0.0, 0, [[], []], [[], []], [[], []]
	else:
		cycles = [[], []] # cycles, paths
		cycle_weights = [[], []] # cycles, paths
		path_constraints_satisfied = [[], []] # cycles, paths
		path_constraints_satisfied_set = set([])

		sol_z = m.getAttr('X', z)
		sol_w = m.getAttr('X', w)
		sol_d = m.getAttr('X', d)
		sol_r = m.getAttr('X', r)		
		total_weights_included = 0.0
		for i in range(k):
			if sol_z[i] >= 0.9:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCycle/Path %d exists; CN = %f." %(i, sol_w[i]))
				sol_x = m.getAttr('X', x)
				sol_c = m.getAttr('X', c)
				cycle_flag = -1
				for ci in range(len(sol_c)):
					if ci % k == i and sol_c[ci] >= 0.9:
						cycle_flag = ci // k
						break
				if cycle_flag == -1:
					cycle = dict()
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if xi % k == i and sol_x[xi] >= 0.9:
							xi_ = xi // k
							x_xi = int(round(sol_x[xi]))
							if xi_ < lseg:
								cycle[('e', xi_)] = x_xi
							elif xi_ < lseg + lc:
								cycle[('c', xi_ - lseg)] = x_xi
							elif xi_ < lseg + lc + ld:
								cycle[('d', xi_ - lseg - lc)] = x_xi
							elif xi_ < lseg + lc + ld + 2 * lsrc:
								assert x_xi == 1
								if (xi_ - lseg - lc - ld) % 2 == 0:
									cycle[('s', (xi_ - lseg - lc - ld) // 2)] = 1 # source edge connected to s
								else:
									cycle[('t', (xi_ - lseg - lc - ld - 1) // 2)] = 1 # source edge connected to t
							else:
								assert x_xi == 1
								if (xi_ - lseg - lc - ld - 2 * lsrc) % 2 == 0:
									nsi = (xi_ - lseg - lc - ld - 2 * lsrc) // 2
									cycle[('ns', nsi)] = 1 # source edge connected to s
								else:
									nti = (xi_ - lseg - lc - ld - 2 * lsrc - 1) // 2
									cycle[('nt', nti)] = 1 # source edge connected to t
					for pi in range(len(pc_list)):
						if sol_r[pi * k + i] >= 0.9:
							path_constraints_s.append(pi)
					cycles[1].append(cycle)
					cycle_weights[1].append(sol_w[i])
					path_constraints_satisfied[1].append(path_constraints_s)
					path_constraints_satisfied_set |= set(path_constraints_s)
				else:
					cycle = dict()
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if xi % k == i and sol_x[xi] >= 0.9:
							xi_ = xi // k
							x_xi = int(round(sol_x[xi]))
							if xi_ < lseg:
								cycle[('e', xi_)] = x_xi
							elif xi_ < lseg + lc:
								cycle[('c', xi_ - lseg)] = x_xi
							elif xi_ < lseg + lc + ld:
								cycle[('d', xi_ - lseg - lc)] = x_xi
							else:
								logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
										"\tError: Cyclic path cannot connect to source nodes.")
								logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tAborted.")
								os.abort()
					for pi in range(len(pc_list)):
						if sol_r[pi * k + i] >= 0.9:
							path_constraints_s.append(pi)
					cycles[0].append(cycle)
					cycle_weights[0].append(sol_w[i])
					path_constraints_satisfied[0].append(path_constraints_s)
					path_constraints_satisfied_set |= set(path_constraints_s)
				for seqi in range(lseg):
					total_weights_included += (sol_x[seqi * k + i] * sol_w[i] * g.sequence_edges[seqi][-2])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
				"Total length weighted CN from cycles/paths = %f/%f." %(total_weights_included, total_weights))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
				"Total num subpath constraints satisfied = %d/%d." %(len(path_constraints_satisfied_set), len(pc_list)))
		return m.Status, total_weights_included, len(path_constraints_satisfied_set), cycles, cycle_weights, path_constraints_satisfied


def minimize_cycles_post(amplicon_id, g, total_weights, node_order, pc_list, init_sol, max_seq_repeat = 2, 
		p_total_weight = 0.9, resolution = 0.1, num_threads = -1, time_limit = 7200, model_prefix = ""):
	"""
	Cycle decomposition by postprocessing the greedy solution

	g: breakpoint graph (object)
	total_weights: float, total length-weighted CN in breakpoint graph g
	node_order: dict maps each node in the input breakpoint graphg to a distinct integer, indicating a total order of the nodes in g
	pc_list: list of subpath constraints to be satisfied, each as a dict that maps an edge to its multiplicity
	init_sol: initial solution returned by maximize_weights_greedy
	max_seq_repeat: integer, maximum multiplicity (num occurrence) allowed for each sequence edge in each cycle/path, 
			default value is 2
	p_total_weight: float between (0, 1), minimum proportion of length-weighted CN to be covered by the resulting cycles or paths, 
			default value is 0.9
	resolution: float, minimum CN for each cycle or path, default value is 0.1
	num_threads: integer, number of working threads for gurobipy, by default it tries to use up all available cores
	time_limit: integer, maximum allowed running time, in seconds, default is 7200 (2 hour)
	model_prefix: output prefix for gurobi *.lp model

	Returns: (1) Status of gurobi optimization model (usually 2 - optimal; 3 - infeasible; 9 - suboptimal/reached time limit)
		(2) Total length weighted CN in resulting cycles/paths
		(3) Total num subpath constraints satisfied by resulting cycles/paths
		(4) List of cycles, each as a dict which maps an edge to its multiplicity in the cycle
		(5) List of the corresponding CN of the above cycles
		(6) Subpath constraints (indices) satisfied by each cycle
		(7) List of paths, each as a dict which maps an edge to its multiplicity in the path
		(8) List of the corresponding CN of the above paths
		(9) Subpath constraints (indices) satisfied by each path
	"""
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Cycle decomposition with initial solution from the greedy strategy.")

	k = len(init_sol[0][0]) + len(init_sol[0][1])
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tReset k (num cycles) to %d." %k)
	p_path_constraints = 0
	path_constraint_indices_ = []
	for paths in (init_sol[2][0] + init_sol[2][1]):
		for pathi in paths:
			if pathi not in path_constraint_indices_:
				path_constraint_indices_.append(pathi)
	if len(pc_list) > 0:
		p_path_constraints = len(path_constraint_indices_) * 0.9999 / len(pc_list)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tRequired proportion of subpath constraints to be satisfied: %f." %p_path_constraints)
	else:
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tProceed without subpath constraints.")

	lseg = len(g.sequence_edges)
	lc = len(g.concordant_edges)
	ld = len(g.discordant_edges)
	lsrc = len(g.source_edges)
	nnodes = len(g.nodes)
	nedges = lseg + lc + ld + 2 * lsrc + 2 * len(g.endnodes)
	endnode_list = [node for node in g.endnodes.keys()]
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tNum nodes to be used in QP = %d." %nnodes)
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tNum edges to be used in QP = %d." %nedges)

	# Gurobi model
	m = gp.Model(model_prefix + "_amplicon" + str(amplicon_id) + "_cycle_postprocessing_k=" + str(k))
		
	# z[i]: indicating whether cycle or path i exists
	z = m.addVars(k, vtype = GRB.BINARY, name = ["z" + str(i) for i in range(k)])
		
	# w[i]: the weight of cycle or path i, continuous variable
	w = m.addVars(k, lb = 0.0, ub = g.max_cn, vtype = GRB.CONTINUOUS, name = ["w" + str(i) for i in range(k)])
		
	# Relationship between w[i] and z[i]
	for i in range(k):
		m.addConstr(w[i] <= z[i] * g.max_cn)
		m.addConstr(w[i] >= z[i] * resolution)

	# x: the number of times an edge occur in cycle or path i
	x_names = []
	for ei in range(nedges):
		for i in range(k):
			x_names.append("x" + str(ei) + "," + str(i))
	x = m.addVars(k * nedges, lb = 0.0, ub = 10.0, vtype = GRB.INTEGER, name = x_names)

	# r and R: subpath constraints
	r = []
	R = []
	if len(pc_list) > 0:
		r_names = []
		R_names = []
		for pi in range(len(pc_list)):
			R_names.append("R" + str(pi))
			for i in range(k):
				r_names.append("r" + str(pi) + "," + str(i))
		r = m.addVars(k * len(pc_list), vtype = GRB.BINARY, name = r_names)
		R = m.addVars(len(pc_list), vtype = GRB.BINARY, name = R_names)

	# Objective: minimize the total number of cycles
	obj = gp.QuadExpr(1.0)
	if len(pc_list) > 0:
		obj = gp.QuadExpr(2.0)
	for i in range(k):
		obj += z[i]
		for seqi in range(lseg):
			obj -= (x[seqi * k + i] * w[i] * g.sequence_edges[seqi][-2] / total_weights)
	for pi in range(len(pc_list)):
		obj -= (R[pi] / len(pc_list))
	m.setObjective(obj, GRB.MINIMIZE)

	# Must include at least 0.9 * total CN weights
	total_weights_expr = gp.QuadExpr(0.0)
	for i in range(k):
		for seqi in range(lseg):
			total_weights_expr += (x[seqi * k + i] * w[i] * g.sequence_edges[seqi][-2])
	m.addConstr(total_weights_expr >= p_total_weight * total_weights)

	# Eulerian constraint
	for node in g.nodes.keys():
		if node in g.endnodes:
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
					ec_expr -= x[(lseg + lc + di) * k + i]
				for srci in g.nodes[node][3]:
					ec_expr -= x[(lseg + lc + ld + 2 * srci) * k + i] # connected to s
					ec_expr -= x[(lseg + lc + ld + 2 * srci) * k + k + i] # connected to t
				m.addConstr(ec_expr == 0.0)
	for i in range(k):
		path_expr = gp.LinExpr(0.0)
		for enodei in range(len(g.endnodes)):
			path_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * enodei) * k + i] # (s, v)
			path_expr -= x[(lseg + lc + ld + 2 * lsrc + 2 * enodei) * k + k + i] # (v, t)
		for srci in range(lsrc): 
			path_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
			path_expr -= x[(lseg + lc + ld + 2 * srci) * k + k + i] # (v, t)
		m.addConstr(path_expr == 0.0)

	# CN constraint
	for seqi in range(lseg):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[seqi * k + i]
		m.addQConstr(cn_expr <= g.sequence_edges[seqi][-1])
	for ci in range(lc):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[(lseg + ci) * k + i]
		m.addQConstr(cn_expr <= g.concordant_edges[ci][-1])
	for di in range(ld):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[(lseg + lc + di) * k + i]
		m.addQConstr(cn_expr <= g.discordant_edges[di][-1])
	for srci in range(lsrc):
		cn_expr = gp.QuadExpr(0.0)
		for i in range(k):
			cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + i]
			cn_expr += w[i] * x[(lseg + lc + ld + 2 * srci) * k + k + i]
		m.addQConstr(cn_expr <= g.source_edges[srci][-1])

	# Occurrence of breakpoints in each cycle/path
	for i in range(k):
		for seqi in range(lseg):
			m.addConstr(x[seqi * k + i] <= max_seq_repeat)
	for i in range(k):
		sum_inv = gp.LinExpr(1.0)
		for di in range(ld):
			dedge = g.discordant_edges[di]
			if dedge[2] == dedge[5]:
				sum_inv += x[(lseg + lc + di) * k + i]
		for di in range(ld):
			m.addConstr(x[(lseg + lc + di) * k + i] <= 2)
			m.addConstr(x[(lseg + lc + di) * k + i] <= sum_inv)

	# c: decomposition i is a cycle, and start at particular node
	c_names = []
	for ni in range(nnodes):
		for i in range(k):
			c_names.append("c" + str(ni) + "," + str(i))
	c = m.addVars(k * nnodes, vtype = GRB.BINARY, name = c_names)
		
	# Relationship between c and x
	for i in range(k):
		cycle_expr = gp.LinExpr(0.0)
		for ni in range(nnodes):
			cycle_expr += c[ni * k + i]
		for eni in range(len(g.endnodes)):
			cycle_expr += x[(lseg + lc + ld + 2 * lsrc + 2 * eni) * k + i] # (s, v)
		for srci in range(lsrc): 
			cycle_expr += x[(lseg + lc + ld + 2 * srci) * k + i] # (s, v)
		m.addConstr(cycle_expr <= 1.0)

	# special request for c added for max_seq_repeat >= 2
	if max_seq_repeat >= 2:
		for i in range(k):
			expr_xc = gp.QuadExpr(0.0)
			for node in g.nodes.keys():
				for ci in set(g.nodes[node][1]):
					expr_xc += c[k * node_order[node] + i] * x[(lseg + ci) * k + i]
				for di in set(g.nodes[node][2]):	
					expr_xc += c[k * node_order[node] + i] * x[(lseg + lc + di) * k + i]
			m.addConstr(expr_xc <= 1.0)

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
			
	# Relationship between c and d
	for i in range(k):
		for ni in range(nnodes):
			m.addConstr(d[ni * k + i] >= c[ni * k + i])
	for i in range(k):
		d_expr = gp.LinExpr(0.0)
		for ni in range(nnodes):
			d_expr += c[ni * k + i]
		d_expr += d[k * nnodes + i] # d_s,i
		m.addConstr(d_expr <= 1.0)
			
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
			dedge = g.discordant_edges[di]
			if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
				m.addConstr(y1[(lseg + lc + di) * k + i] == 0)
				m.addConstr(y2[(lseg + lc + di) * k + i] == 0)
	for i in range(k):
		t_expr_x = gp.LinExpr(0.0)
		t_expr_y = gp.LinExpr(0.0)
		t_expr_yd = gp.QuadExpr(0.0)
		for node in endnode_list:
			t_expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i]
			t_expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i]
			t_expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i] * \
					(d[k * nnodes + k + i] - d[k * node_order[node] + i]) # node -> t
			expr_x = gp.LinExpr(0.0)
			expr_y = gp.LinExpr(0.0)
			expr_xc = gp.QuadExpr(0.0)
			expr_yd = gp.QuadExpr(0.0)
			for seqi in g.nodes[node][0]:
				sseg = g.sequence_edges[seqi]
				node_ = (sseg[0], sseg[1], '-')
				if node_ == node:
					node_ = (sseg[0], sseg[2], '+')
				expr_x += x[seqi * k + i]
				expr_xc += x[seqi * k + i] * c[k * node_order[node] + i]
				if node_order[node_] <= node_order[node]:
					expr_y += y1[seqi * k + i]
					expr_yd += y1[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				else:
					expr_y += y2[seqi * k + i]
					expr_yd += y2[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
						
			expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + i] # from s
			expr_xc += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + i] * c[k * node_order[node] + i]
			expr_x += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i] # to t
			expr_xc += x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + k + i] * c[k * node_order[node] + i]
			expr_y += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + i] # from s
			expr_yd += y1[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) * k + i] * \
					(d[k * node_order[node] + i] - d[k * nnodes + i])
			m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
			m.addConstr(expr_y <= 1.0)
			m.addConstr(expr_y * nedges * k + expr_xc >= expr_x)
			m.addConstr(expr_yd * nedges * k + expr_xc >= expr_x)
			
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
				expr_xc = gp.QuadExpr(0.0)
				expr_yd = gp.QuadExpr(0.0)
				for seqi in g.nodes[node][0]:
					sseg = g.sequence_edges[seqi]
					node_ = (sseg[0], sseg[1], '-')
					if node_ == node:
						node_ = (sseg[0], sseg[2], '+')
					expr_x += x[seqi * k + i]
					expr_xc += x[seqi * k + i] * c[k * node_order[node] + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[seqi * k + i]
						expr_yd += y1[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[seqi * k + i]
						expr_yd += y2[seqi * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				for ci in g.nodes[node][1]:
					cedge = g.concordant_edges[ci]
					node_ = (cedge[0], cedge[1], cedge[2])
					if node_ == node:
						node_ = (cedge[3], cedge[4], cedge[5])
					expr_x += x[(lseg + ci) * k + i]
					expr_xc += x[(lseg + ci) * k + i] * c[k * node_order[node] + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[(lseg + ci) * k + i]
						expr_yd += y1[(lseg + ci) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[(lseg + ci) * k + i]
						expr_yd += y2[(lseg + ci) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				for di in g.nodes[node][2]:
					dedge = g.discordant_edges[di]
					node_ = (dedge[0], dedge[1], dedge[2])
					if node_ == node:
						node_ = (dedge[3], dedge[4], dedge[5])
					expr_x += x[(lseg + lc + di) * k + i]
					expr_xc += x[(lseg + lc + di) * k + i] * c[k * node_order[node] + i]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[(lseg + lc + di) * k + i]
						expr_yd += y1[(lseg + lc + di) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
					else:
						expr_y += y2[(lseg + lc + di) * k + i]
						expr_yd += y2[(lseg + lc + di) * k + i] * (d[k * node_order[node] + i] - d[k * node_order[node_] + i])
				for srci in g.nodes[node][3]:
					expr_x += x[(lseg + lc + ld + 2 * srci) * k + i]
					expr_x += x[(lseg + lc + ld + 2 * srci) * k + k + i]
					expr_xc += x[(lseg + lc + ld + 2 * srci) * k + i] * c[k * node_order[node] + i]
					expr_xc += x[(lseg + lc + ld + 2 * srci) * k + k + i] * c[k * node_order[node] + i]
					expr_y += y1[(lseg + lc + ld + 2 * srci) * k + i]
					expr_yd += y1[(lseg + lc + ld + 2 * srci) * k + i] * \
							(d[k * node_order[node] + i] - d[k * nnodes + i])
				m.addConstr(expr_x * (nnodes + 2) >= d[k * node_order[node] + i])
				m.addConstr(expr_y <= 1.0)
				m.addConstr(expr_y * nedges * k + expr_xc >= expr_x)
				m.addConstr(expr_yd * nedges * k + expr_xc >= expr_x)

	# Subpath constraints
	if len(pc_list) > 0:
		sum_R = gp.LinExpr(0.0)
		for pi in range(len(pc_list)):
			sum_R += R[pi]
			path_constraint_ = pc_list[pi]
			sum_r = gp.LinExpr(0.0)
			for ri in range(pi * k, (pi + 1) * k):
				sum_r += r[ri]
				m.addConstr(R[pi] >= r[ri])
			m.addConstr(sum_r >= R[pi])	
			for edge in path_constraint_.keys():
				for i in range(k):
					if edge[0] == 's':
						m.addConstr(x[edge[1] * k + i] >= r[pi * k + i] * path_constraint_[edge])
					elif edge[0] == 'c':
						m.addConstr(x[(lseg + edge[1]) * k + i] >= r[pi * k + i] * path_constraint_[edge])
					else:
						m.addConstr(x[(lseg + lc + edge[1]) * k + i] >= r[pi * k + i] * path_constraint_[edge])
		m.addConstr(sum_R >= p_path_constraints * len(pc_list))			

	# Initialize variables
	for i in range(len(init_sol[0][0])):
		z[i].start = 1
		w[i].start = init_sol[1][0][i]
		for (v, vi) in init_sol[0][0][i].keys():
			if v == 'x':
				x[vi * k + i].start = init_sol[0][0][i][(v, vi)]
			elif v == 'c':
				c[vi * k + i].start = init_sol[0][0][i][(v, vi)]
			elif v == 'd':
				d[vi * k + i].start = init_sol[0][0][i][(v, vi)]
			elif v == 'y1':
				y1[vi * k + i].start = init_sol[0][0][i][(v, vi)]
			elif v == 'y2':
				y2[vi * k + i].start = init_sol[0][0][i][(v, vi)]
	for i in range(len(init_sol[0][1])):
		i_ = i + len(init_sol[0][0])
		z[i_].start = 1
		w[i_].start = init_sol[1][1][i]
		for (v, vi) in init_sol[0][1][i].keys():
			if v == 'x':
				x[vi * k + i_].start = init_sol[0][1][i][(v, vi)]
			elif v == 'c':
				c[vi * k + i_].start = init_sol[0][1][i][(v, vi)]
			elif v == 'd':
				d[vi * k + i_].start = init_sol[0][1][i][(v, vi)]
			elif v == 'y1':
				y1[vi * k + i_].start = init_sol[0][1][i][(v, vi)]
			elif v == 'y2':
				y2[vi * k + i_].start = init_sol[0][1][i][(v, vi)]
	for i in range(len(init_sol[2][0])):
		for pi in init_sol[2][0][i]:
			r[pi * k + i].start = 1
			R[pi].start = 1
	for i in range(len(init_sol[2][1])):
		i_ = i + len(init_sol[2][0])
		for pi in init_sol[2][1][i]:
			r[pi * k + i_].start = 1
			R[pi].start = 1
	m.update()

	m.setParam(GRB.Param.LogToConsole, 0)
	if num_threads > 0:
		m.setParam(GRB.Param.Threads, num_threads)
	m.setParam(GRB.Param.NonConvex, 2)
	m.setParam(GRB.Param.Heuristics, 0.25)
	m.setParam(GRB.Param.TimeLimit, max(time_limit, ld * 300)) # each breakpoint edge is assigned 5 minutes
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCompleted gurobi model setup.")
	lp_fn = model_prefix + "_amplicon" + str(amplicon_id) + "_postprocessing_model.lp"
	m.write(lp_fn)
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tWrote model to file: %s." %lp_fn)
	log_fn = lp_fn[:-2] + "log"
	m.setParam(GRB.Param.LogFile, log_fn)
	m.optimize()
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCompleted optimization with status %d." %m.Status)
	
	if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tModel is infeasible.")
		return GRB.INFEASIBLE, 0.0, 0, [[], []], [[], []], [[], []]
	else:
		cycles = [[], []] # cycles, paths
		cycle_weights = [[], []] # cycles, paths
		path_constraints_satisfied = [[], []] # cycles, paths
		path_constraints_satisfied_set = set([])

		sol_z = m.getAttr('X', z)
		sol_w = m.getAttr('X', w)
		sol_d = m.getAttr('X', d)
		sol_r = m.getAttr('X', r)
		total_weights_included = 0.0
		for i in range(k):
			if sol_z[i] >= 0.9:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCycle/Path %d exists; CN = %f." %(i, sol_w[i]))
				sol_x = m.getAttr('X', x)
				sol_c = m.getAttr('X', c)
				cycle_flag = -1
				for ci in range(len(sol_c)):
					if ci % k == i and sol_c[ci] >= 0.9:
						cycle_flag = ci // k
						break
				if cycle_flag == -1:
					cycle = dict()
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if xi % k == i and sol_x[xi] >= 0.9:
							xi_ = xi // k
							x_xi = int(round(sol_x[xi]))
							if xi_ < lseg:
								cycle[('e', xi_)] = x_xi
							elif xi_ < lseg + lc:
								cycle[('c', xi_ - lseg)] = x_xi
							elif xi_ < lseg + lc + ld:
								cycle[('d', xi_ - lseg - lc)] = x_xi
							elif xi_ < lseg + lc + ld + 2 * lsrc:
								assert x_xi == 1
								if (xi_ - lseg - lc - ld) % 2 == 0:
									cycle[('s', (xi_ - lseg - lc - ld) // 2)] = 1 # source edge connected to s
								else:
									cycle[('t', (xi_ - lseg - lc - ld - 1) // 2)] = 1 # source edge connected to t
							else:
								assert x_xi == 1
								if (xi_ - lseg - lc - ld - 2 * lsrc) % 2 == 0:
									nsi = (xi_ - lseg - lc - ld - 2 * lsrc) // 2
									cycle[('ns', nsi)] = 1 # source edge connected to s
								else:
									nti = (xi_ - lseg - lc - ld - 2 * lsrc - 1) // 2
									cycle[('nt', nti)] = 1 # source edge connected to t
					for pi in range(len(pc_list)):
						if sol_r[pi * k + i] >= 0.9:
							path_constraints_s.append(pi)
					cycles[1].append(cycle)
					cycle_weights[1].append(sol_w[i])
					path_constraints_satisfied[1].append(path_constraints_s)
					path_constraints_satisfied_set |= set(path_constraints_s)
				else:
					cycle = dict()
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if xi % k == i and sol_x[xi] >= 0.9:
							xi_ = xi // k
							x_xi = int(round(sol_x[xi]))
							if xi_ < lseg:
								cycle[('e', xi_)] = x_xi
							elif xi_ < lseg + lc:
								cycle[('c', xi_ - lseg)] = x_xi
							elif xi_ < lseg + lc + ld:
								cycle[('d', xi_ - lseg - lc)] = x_xi
							else:
								logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
										"\tError: Cyclic path cannot connect to source nodes.")
								logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tAborted.")
								os.abort()
					for pi in range(len(pc_list)):
						if sol_r[pi * k + i] >= 0.9:
							path_constraints_s.append(pi)
					cycles[0].append(cycle)
					cycle_weights[0].append(sol_w[i])
					path_constraints_satisfied[0].append(path_constraints_s)
					path_constraints_satisfied_set |= set(path_constraints_s)
				for seqi in range(lseg):
					total_weights_included += (sol_x[seqi * k + i] * sol_w[i] * g.sequence_edges[seqi][-2])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
				"Total length weighted CN from cycles/paths = %f/%f." %(total_weights_included, total_weights))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
				"Total num subpath constraints satisfied = %d/%d." %(len(path_constraints_satisfied_set), len(pc_list)))
		return m.Status, total_weights_included, len(path_constraints_satisfied_set), cycles, cycle_weights, path_constraints_satisfied


def maximize_weights_greedy(amplicon_id, g, total_weights, node_order, pc_list, alpha = 0.01, max_seq_repeat = 2,
			p_total_weight = 0.9, resolution = 0.1, cn_tol = 0.005, p_subpaths = 0.9, num_threads = -1, 
			postprocess = 0, time_limit = 7200, model_prefix = ""):
	"""
	Greedy cycle decomposition by maximizing the length-weighted CN of a single cycle/path

	amplicon_id: integer, amplicon ID
	g: breakpoint graph (object)
	total_weights: float, total length-weighted CN in breakpoint graph g
	node_order: dict maps each node in the input breakpoint graphg to a distinct integer, indicating a total order of the nodes in g
	pc_list: list of subpath constraints to be satisfied, each as a dict that maps an edge to its multiplicity
	alpha: float, parameter for multi-objective optimization, default value is 0.01
		maximizing total length-weighted CN + 
		(alpha * remaining length-weighted CN in the graph / num remaining unsatisfied subpath constraints) * 
		num subpath constraints satisfied by the next cycle or path
		*** when alpha < 0, just maximizing total length-weighted CN
	max_seq_repeat: integer, maximum multiplicity (num occurrence) allowed for each sequence edge in each cycle/path, 
			default value is 2
	p_total_weight: float between (0, 1), minimum proportion of length-weighted CN to be covered by the resulting cycles or paths, 
			default value is 0.9
	resolution: float, minimum CN for each cycle or path, default value is 0.1
	cn_tol: float between (0, 1), terminate greedy cycle decomposition when the next cycle/path has total length weighted CN 
		< cn_tol * total_weights, default value is 0.005
	p_subpaths: float between (0, 1), minimum proportion of subpath constraints to be satisfied by the resulting cycles or paths, 
		default value is 0.9
	num_threads: integer, number of working threads for gurobipy, by default it tries to use up all available cores
	time_limit: integer, maximum allowed running time, in seconds, default is 7200 (2 hour)
	model_prefix: output prefix for gurobi *.lp model

	Returns: (1) Total length weighted CN in resulting cycles/paths
		(2) Total num subpath constraints satisfied by resulting cycles/paths
		(3) List of cycles, each as a dict which maps an edge to its multiplicity in the cycle
		(4) List of the corresponding CN of the above cycles
		(5) Subpath constraints (indices) satisfied by each cycle
		(6) List of paths, each as a dict which maps an edge to its multiplicity in the path
		(7) List of the corresponding CN of the above paths
		(8) Subpath constraints (indices) satisfied by each path 
	"""
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Integer program too large, perform greedy cycle decomposition.")
	lseg = len(g.sequence_edges)
	lc = len(g.concordant_edges)
	ld = len(g.discordant_edges)
	lsrc = len(g.source_edges)
	nnodes = len(g.nodes)
	nedges = lseg + lc + ld + 2 * lsrc + 2 * len(g.endnodes)
	endnode_list = [node for node in g.endnodes.keys()]
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tNum nodes to be used in QP = %d." %nnodes)
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tNum edges to be used in QP = %d." %nedges)
		
	remaining_weights = total_weights
	unsatisfied_pc = [i for i in range(len(pc_list))]
	remaining_CN = dict()
	for segi in range(lseg):
		remaining_CN[('s', segi)] = g.sequence_edges[segi][-1]
	for ci in range(lc):
		remaining_CN[('c', ci)] = g.concordant_edges[ci][-1]
	for di in range(ld):
		remaining_CN[('d', di)] = g.discordant_edges[di][-1]
	for srci in range(lsrc):
		remaining_CN[('src', srci)] = g.source_edges[srci][-1]
	next_w = resolution * 1.1
	cycle_id = 0
	num_unsatisfied_pc = len(pc_list)
	cycles = [[], []] # cycles, paths
	cycle_weights = [[], []] # cycles, paths
	path_constraints_satisfied = [[], []] # cycles, paths
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
			"\tGreedy cycle decomposition with length weighted CN = %f and num subpath constraints = %d." 
			%(remaining_weights, num_unsatisfied_pc))	

	while next_w >= resolution and (remaining_weights > (1.0 - p_total_weight) * total_weights or \
		num_unsatisfied_pc > math.floor((1.0 - p_subpaths) * len(pc_list))):
		pp = 1.0
		if alpha > 0 and num_unsatisfied_pc > 0:
			pp = alpha * remaining_weights / num_unsatisfied_pc # multi - objective optimization parameter
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tIteration %d." %(cycle_id + 1))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
				"\tRemaining length weighted CN: %f/%f." %(remaining_weights, total_weights))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
				"\tRemaining subpath constraints to be satisfied: %d/%d." %(num_unsatisfied_pc, len(pc_list)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tMultiplication factor for subpath constraints: %f." %pp)

		# Gurobi model
		m = gp.Model(model_prefix + "_amplicon" + str(amplicon_id) + "_cycle_decomposition_greedy_" + str(cycle_id + 1))

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
			if node in g.endnodes:
				m.addConstr(x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node))] + \
					x[(lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)) + 1] == x[g.nodes[node][0][0]])
			else:
				ec_expr = gp.LinExpr(0.0)
				for seqi in g.nodes[node][0]:
					ec_expr += x[seqi]
				for ci in g.nodes[node][1]:
					ec_expr -= x[lseg + ci]
				for di in g.nodes[node][2]:
					ec_expr -= x[lseg + lc + di]
				for srci in g.nodes[node][3]:
					ec_expr -= x[lseg + lc + ld + 2 * srci] # connected to s
					ec_expr -= x[lseg + lc + ld + 2 * srci + 1] # connected to t
				m.addConstr(ec_expr == 0.0)
		path_expr = gp.LinExpr(0.0)
		for enodei in range(len(g.endnodes)):
			path_expr += x[lseg + lc + ld + 2 * lsrc + 2 * enodei] # (s, v)
			path_expr -= x[lseg + lc + ld + 2 * lsrc + 2 * enodei + 1] # (v, t)
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
			if g.discordant_edges[di][-1] < resolution:
				m.addConstr(x[lseg + lc + di] == 0.0)
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
						"\tIgnored discordant edge at index %d." %(di))
		for srci in range(lsrc):
			cn_expr = gp.QuadExpr(0.0)
			cn_expr += w[0] * x[lseg + lc + ld + 2 * srci]
			cn_expr += w[0] * x[lseg + lc + ld + 2 * srci + 1]
			m.addQConstr(cn_expr <= remaining_CN[('src', srci)])
			
		# Occurrence of breakpoints in each cycle/path
		for seqi in range(lseg):
			m.addConstr(x[seqi] <= max_seq_repeat)
		sum_inv = gp.LinExpr(1.0)
		for di in range(ld):
			dedge = g.discordant_edges[di]
			if dedge[2] == dedge[5]:
				sum_inv += x[lseg + lc + di]
		for di in range(ld):
			m.addConstr(x[lseg + lc + di] <= 2)
			m.addConstr(x[lseg + lc + di] <= sum_inv)
			
		# c: decomposition i is a cycle, and start at particular node
		c_names = []
		for ni in range(nnodes):
			c_names.append("c" + str(ni))
		c = m.addVars(nnodes, vtype = GRB.BINARY, name = c_names)

		# Relationship between c and x
		cycle_expr = gp.LinExpr(0.0)
		for ni in range(nnodes):
			cycle_expr += c[ni]
		for eni in range(len(g.endnodes)):
			cycle_expr += x[lseg + lc + ld + 2 * lsrc + 2 * eni] # (s, v)
		for srci in range(lsrc): 
			cycle_expr += x[lseg + lc + ld + 2 * srci] # (s, v)
		m.addConstr(cycle_expr <= 1.0)

		# special request for c added for max_seq_repeat >= 2
		if max_seq_repeat >= 2:
			expr_xc = gp.QuadExpr(0.0)
			for node in g.nodes.keys():
				for ci in set(g.nodes[node][1]):
					expr_xc += c[node_order[node]] * x[lseg + ci]
				for di in set(g.nodes[node][2]):	
					expr_xc += c[node_order[node]] * x[lseg + lc + di]
			m.addConstr(expr_xc <= 1.0)

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
			dedge = g.discordant_edges[di]
			if dedge[0] == dedge[3] and dedge[1] == dedge[4] and dedge[2] == dedge[5]: # exclude self loops
				m.addConstr(y1[lseg + lc + di] == 0)
				m.addConstr(y2[lseg + lc + di] == 0)
		t_expr_x = gp.LinExpr(0.0)
		t_expr_y = gp.LinExpr(0.0)
		t_expr_yd = gp.QuadExpr(0.0)
		for node in endnode_list:
			t_expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node) + 1]
			t_expr_y += y1[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node) + 1]
			t_expr_yd += y1[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node) + 1] * \
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
						
			expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)] # from s
			expr_xc += x[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)] * c[node_order[node]]
			expr_x += x[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node) + 1] # to t
			expr_xc += x[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node) + 1] * c[node_order[node]]
			expr_y += y1[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)] # from s
			expr_yd += y1[lseg + lc + ld + 2 * lsrc + 2 * endnode_list.index(node)] * \
								(d[node_order[node]] - d[nnodes])
			m.addConstr(expr_x * (nnodes + 2) >= d[node_order[node]])
			m.addConstr(expr_y <= 1.0)
			m.addConstr(expr_y * nedges + expr_xc >= expr_x)
			m.addConstr(expr_yd * nedges + expr_xc >= expr_x)
		for srci in range(lsrc):
			srce = g.source_edges[srci]
			t_expr_x += x[lseg + lc + ld + 2 * srci + 1]
			t_expr_y += y1[lseg + lc + ld + 2 * srci + 1]
			t_expr_yd += y1[lseg + lc + ld + 2 * srci + 1] * (d[nnodes + 1] - d[node_order[(srce[3], srce[4], srce[5])]])
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
					cedge = g.concordant_edges[ci]
					node_ = (cedge[0], cedge[1], cedge[2])
					if node_ == node:
						node_ = (cedge[3], cedge[4], cedge[5])
					expr_x += x[lseg + ci]
					expr_xc += x[lseg + ci] * c[node_order[node]]
					if node_order[node_] <= node_order[node]:
						expr_y += y1[lseg + ci]
						expr_yd += y1[lseg + ci] * (d[node_order[node]] - d[node_order[node_]])
					else:
						expr_y += y2[lseg + ci]
						expr_yd += y2[lseg + ci] * (d[node_order[node]] - d[node_order[node_]])
				for di in g.nodes[node][2]:
					dedge = g.discordant_edges[di]
					node_ = (dedge[0], dedge[1], dedge[2])
					if node_ == node:
						node_ = (dedge[3], dedge[4], dedge[5])
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
		
		m.setParam(GRB.Param.LogToConsole, 0)
		if num_threads > 0:
			m.setParam(GRB.Param.Threads, num_threads)
		m.setParam(GRB.Param.NonConvex, 2)
		m.setParam(GRB.Param.TimeLimit, time_limit) # each breakpoint edge is assigned 5 minutes
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCompleted gurobi setup for model %d." %(cycle_id + 1))
		lp_fn = model_prefix + "_amplicon" + str(amplicon_id) + "_greedy_model_" + str(cycle_id + 1) + "_alpha=" + str(alpha) + ".lp"
		m.write(lp_fn)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tWrote model to %d file: %s." %(cycle_id + 1, lp_fn))
		log_fn = lp_fn[:-2] + "log"
		m.setParam(GRB.Param.LogFile, log_fn)
		m.optimize()
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCompleted optimization of model %d with status %d." %(cycle_id + 1, m.Status))
		
		if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tModel %d is infeasible." %(cycle_id + 1))
			logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tIteration terminated.")
			break
		else:
			cycle_id += 1
			sol_z = m.getAttr('X', z)
			sol_w = m.getAttr('X', w)
			sol_d = m.getAttr('X', d)
			sol_r = m.getAttr('X', r)
			total_weights_included = 0.0
			if sol_z[0] >= 0.9:
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tNext cycle/path exists; CN = %f." %(sol_w[0]))
				next_w = sol_w[0]
				if next_w < resolution:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + 
							"\tCN less than resolution, iteration terminated successfully.")
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
				cycle_for_postprocess = dict()
				for xi in range(len(sol_x)):
					cycle_for_postprocess[('x', xi)] = sol_x[xi]
				for ci in range(len(sol_c)):
					cycle_for_postprocess[('c', ci)] = sol_c[ci]
				for di in range(len(sol_d)):
					cycle_for_postprocess[('d', di)] = sol_d[di]
				for yi in range(len(sol_y1)):
					cycle_for_postprocess[('y1', yi)] = sol_y1[yi]
				for yi in range(len(sol_y2)):
					cycle_for_postprocess[('y2', yi)] = sol_y2[yi]
				if cycle_flag == -1:
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if sol_x[xi] >= 0.9:
							x_xi = int(round(sol_x[xi]))
							if xi < lseg:
								cycle[('e', xi)] = x_xi
								remaining_CN[('s', xi)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('s', xi)] < resolution:
									remaining_CN[('s', xi)] = 0.0
							elif xi < lseg + lc:
								cycle[('c', xi - lseg)] = x_xi
								remaining_CN[('c', xi - lseg)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('c', xi - lseg)] < resolution:
									remaining_CN[('c', xi - lseg)] = 0.0
							elif xi < lseg + lc + ld:
								cycle[('d', xi - lseg - lc)] = x_xi
								remaining_CN[('d', xi - lseg - lc)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('d', xi - lseg - lc)] < resolution:
									remaining_CN[('d', xi - lseg - lc)] = 0.0
							elif xi < lseg + lc + ld + 2 * lsrc:
								assert x_xi == 1
								if (xi - lseg - lc - ld) % 2 == 0:
									cycle[('s', (xi_ - lseg - lc - ld) // 2)] = 1 # source edge connected to s
									remaining_CN[('src', (xi - lseg - lc - ld) // 2)] -= sol_w[0]
									if remaining_CN[('src', (xi - lseg - lc - ld) // 2)] < resolution:
										remaining_CN[('src', (xi - lseg - lc - ld) // 2)] = 0.0
								else:
									cycle[('t', (xi_ - lseg - lc - ld - 1) // 2)] = 1 # source edge connected to t
									remaining_CN[('src', (xi - lseg - lc - ld - 1) // 2)] -= sol_w[0]
									if remaining_CN[('src', (xi - lseg - lc - ld - 1) // 2)] < resolution:
										remaining_CN[('src', (xi - lseg - lc - ld - 1) // 2)] = 0.0
							else:
								assert x_xi == 1
								if (xi_ - lseg - lc - ld - 2 * lsrc) % 2 == 0:
									nsi = (xi_ - lseg - lc - ld - 2 * lsrc) // 2
									cycle[('ns', nsi)] = 1 # source edge connected to s
								else:
									nti = (xi_ - lseg - lc - ld - 2 * lsrc - 1) // 2
									cycle[('nt', nti)] = 1 # source edge connected to t
					for pi in range(len(pc_list)):
						if sol_r[pi] >= 0.9:
							path_constraints_s.append(pi)
							unsatisfied_pc[pi] = -1
					if postprocess == 1:
						cycles[1].append(cycle_for_postprocess)
					else:
						cycles[1].append(cycle)
					cycle_weights[1].append(sol_w[0])
					path_constraints_satisfied[1].append(path_constraints_s)
				else:
					path_constraints_s = []
					for xi in range(len(sol_x)):
						if sol_x[xi] >= 0.9:
							x_xi = int(round(sol_x[xi]))
							if xi < lseg:
								cycle[('e', xi)] = x_xi
								remaining_CN[('s', xi)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('s', xi)] < resolution:
									remaining_CN[('s', xi)] = 0.0
							elif xi < lseg + lc:
								cycle[('c', xi - lseg)] = x_xi
								remaining_CN[('c', xi - lseg)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('c', xi - lseg)] < resolution:
									remaining_CN[('c', xi - lseg)] = 0.0
							elif xi < lseg + lc + ld:
								cycle[('d', xi - lseg - lc)] = x_xi
								remaining_CN[('d', xi - lseg - lc)] -= (sol_x[xi] * sol_w[0])
								if remaining_CN[('d', xi - lseg - lc)] < resolution:
									remaining_CN[('d', xi - lseg - lc)] = 0.0
							else:
								logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
										"\tError: Cyclic path cannot connect to source nodes.")
								logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tAborted.")
								os.abort()
					for pi in range(len(pc_list)):
						if sol_r[pi] >= 0.9:
							path_constraints_s.append(pi)
							unsatisfied_pc[pi] = -1
					if postprocess == 1:
						cycles[1].append(cycle_for_postprocess)
					else:
						cycles[0].append(cycle)
					cycle_weights[0].append(sol_w[0])
					path_constraints_satisfied[0].append(path_constraints_s)
				for seqi in range(lseg):
					total_weights_included += (sol_x[seqi] * sol_w[0] * g.sequence_edges[seqi][-2])
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
					"Total length weighted CN from cycle/path %d = %f/%f." %(cycle_id, total_weights_included, total_weights))
				remaining_weights -= total_weights_included
				if total_weights_included < cn_tol * total_weights:
					num_unsatisfied_pc = 0
					for i in range(len(pc_list)):
						if unsatisfied_pc[i] >= 0:
							num_unsatisfied_pc += 1
					logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + 
						"\tProportion of length weighted CN less than cn_tol, iteration terminated.")
					break
			else:
				break
		num_unsatisfied_pc = 0
		for i in range(len(pc_list)):
			if unsatisfied_pc[i] >= 0:
				num_unsatisfied_pc += 1
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
			"\tRemaining subpath constraints to be satisfied: %d/%d." %(num_unsatisfied_pc, len(pc_list)))
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tProceed to next iteration.")
	if next_w < resolution:
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCycle/path reaches CN resolution, greedy iteration completed.")
	elif (num_unsatisfied_pc <= math.floor((1.0 - p_subpaths) * len(pc_list))) and (remaining_weights <= (1.0 - p_total_weight) * total_weights):
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCycles/paths explain sufficient CN in graph, greedy iteration completed.")
	else:
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + 
				"\tProportion of length weighted CN less than cn_tol, greedy iteration completed.") 
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
			"Total length weighted CN from cycles/paths = %f/%f." %(total_weights - remaining_weights, total_weights))
	logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + \
			"Total num subpath constraints satisfied = %d/%d." %(len(pc_list) - num_unsatisfied_pc, len(pc_list)))
	return total_weights - remaining_weights, len(pc_list) - num_unsatisfied_pc, cycles, cycle_weights, path_constraints_satisfied


def eulerian_cycle_t(g, edges_next_cycle, path_constraints_next_cycle, path_constraints_support):
	"""
	Return an eulerian traversal of a cycle, represented by a dict of edges

	g: breakpoint graph (object)
	edges_next_cycle: subgraph induced by the cycle, as a dict that maps an edge to its multiplicity
	path_constraints_next_cycle: list of subpath constraints to be satisfied, 
		each as a list of alternating nodes and edges
		***
		Note: the traversal may not satisfy all subpath constraints
		in case not all subpath constraints are satisfied, return the eulerian traversal satisfying the
		maximum number of subpath constraints
		***
	path_constraints_support: num long reads supporting each subpath constraint
	"""
	lseg = len(g.sequence_edges)

	eulerian_cycle = [] # A cycle is edge - node list starting and ending with the same edge
				# Since Eulerian, there could be subcycles in the middle of a cycle
	eulerian_cycle_ = [] # Cycle in AA cycle format
	best_cycle = [] # Cycle in AA cycle format
	valid = 0
	num_trials = 0
	l = len(path_constraints_next_cycle)
	unsatisfied_path_metric = [range(l), 100 * l, 100 * max(path_constraints_support + [0])]
	while valid <= 0 and num_trials < 1000:
		valid = 1
		num_trials += 1
		eulerian_cycle = []
		eulerian_cycle_ = []
		edges_cur = edges_next_cycle.copy()
		last_seq_edge = lseg # Start with the edge with smallest index and on the positive strand
		for edge in edges_cur.keys():
			if edge[0] == 'e':
				last_seq_edge = min(last_seq_edge, edge[1])
		last_edge_dir = '+'
		eulerian_cycle.append(('s', last_seq_edge))
		eulerian_cycle_.append(str(last_seq_edge + 1) + '+')
		while len(edges_cur) > 0:
			seq_edge = g.sequence_edges[last_seq_edge]
			node = (seq_edge[0], seq_edge[2], '+')
			if last_edge_dir == '-':
				node = (seq_edge[0], seq_edge[1], '-')
			eulerian_cycle.append(node)
			next_bp_edges = [] # Since cycle, only consider discordant edges and concordant edges
			for ci in g.nodes[node][1]:
				next_bp_edges.append(('c', ci))
			for di in g.nodes[node][2]:
				next_bp_edges.append(('d', di))
			del_list = [i for i in range(len(next_bp_edges)) if next_bp_edges[i] not in edges_cur]
			for i in del_list[::-1]:
				del next_bp_edges[i]
			if len(next_bp_edges) == 0:
				valid = 0
				break
			if len(next_bp_edges) == 1: # No branching on the path
				eulerian_cycle.append(next_bp_edges[0])
				edges_cur[next_bp_edges[0]] = int(edges_cur[next_bp_edges[0]]) - 1
				if edges_cur[next_bp_edges[0]] == 0:
					del edges_cur[next_bp_edges[0]]
				bp_edge = []
				if next_bp_edges[0][0] == 'c':
					bp_edge = g.concordant_edges[next_bp_edges[0][1]][:6]
				else:
					bp_edge = g.discordant_edges[next_bp_edges[0][1]][:6]
				node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
				if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
					node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
				eulerian_cycle.append(node_)
				last_seq_edge = g.nodes[node_][0][0]
				eulerian_cycle.append(('s', last_seq_edge))
				if node_[2] == '-':
					last_edge_dir = '+'
					eulerian_cycle_.append(str(last_seq_edge + 1) + '+')
				else:
					last_edge_dir = '-'
					eulerian_cycle_.append(str(last_seq_edge + 1) + '-')
				edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
				if edges_cur[('e', last_seq_edge)] == 0:
					del edges_cur[('e', last_seq_edge)]
			else:
				r = random.randint(0, len(next_bp_edges) - 1)
				eulerian_cycle.append(next_bp_edges[r])
				edges_cur[next_bp_edges[r]] = int(edges_cur[next_bp_edges[r]]) - 1
				if edges_cur[next_bp_edges[r]] == 0:
					del edges_cur[next_bp_edges[r]]
				bp_edge = []
				if next_bp_edges[r][0] == 'c':
					bp_edge = g.concordant_edges[next_bp_edges[r][1]][:6]
				else:
					bp_edge = g.discordant_edges[next_bp_edges[r][1]][:6]
				node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
				if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
					node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
				eulerian_cycle.append(node_)
				last_seq_edge = g.nodes[node_][0][0]
				eulerian_cycle.append(('s', last_seq_edge))
				if node_[2] == '-':
					last_edge_dir = '+'
					eulerian_cycle_.append(str(last_seq_edge + 1) + '+')
				else:
					last_edge_dir = '-'
					eulerian_cycle_.append(str(last_seq_edge + 1) + '-')
				edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
				if edges_cur[('e', last_seq_edge)] == 0:
					del edges_cur[('e', last_seq_edge)]
		if valid == 1 and len(best_cycle) == 0:
			best_cycle = eulerian_cycle_
		path_metric = [[], 0, 0]
		# check if the remaining path constraints are satisfied
		for pathi in range(len(path_constraints_next_cycle)):
			path_ = path_constraints_next_cycle[pathi]
			path0 = path_[0]
			s = 0
			for ei in range(len(eulerian_cycle) - 1):
				obj = eulerian_cycle[ei]
				if obj == path0:
					s_ = 1
					for i in range(len(path_)):
						if eulerian_cycle[:-1][(ei + i) % (len(eulerian_cycle) - 1)] != path_[i]:
							s_ = 0 
							break
					if s_ == 1:
						s = 1
						break
					else:
						s_ = 1
						for i in range(len(path_)):
							if eulerian_cycle[:-1][ei - i] != path_[i]:
								s_ = 0 
								break
						if s_ == 1:
							s = 1
							break
			if s == 0 and valid == 1:
				path_metric[0].append(pathi)
				path_metric[1] += len(path_)
				path_metric[2] += path_constraints_support[pathi]
		if valid == 1 and len(path_metric[0]) > 0:
			valid = -1
		if valid != 0 and (len(path_metric[0]) < len(unsatisfied_path_metric[0])) or \
			(len(path_metric[0]) == len(unsatisfied_path_metric[0]) and path_metric[1] < unsatisfied_path_metric[1]) or \
			(len(path_metric[0]) == len(unsatisfied_path_metric[0]) and path_metric[1] == unsatisfied_path_metric[1] and \
			path_metric[2] < unsatisfied_path_metric[2]):
			unsatisfied_path_metric[0] = path_metric[0]
			unsatisfied_path_metric[1] = path_metric[1]
			unsatisfied_path_metric[1] = path_metric[1]
			best_cycle = eulerian_cycle_
	if len(unsatisfied_path_metric[0]) == 0:
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tCycle satisfies all subpath constraints.")
	else:
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tThe following path constraints are not satisfied:")
		for pathi in unsatisfied_path_metric[0]:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\t%s" %path_constraints_next_cycle[pathi])
	return best_cycle


def eulerian_path_t(g, edges_next_path, path_constraints_next_path, path_constraints_support):
	"""
	Return an eulerian traversal of an s-t walk, represented by a dict of edges

	g: breakpoint graph (object)
	edges_next_path: subgraph induced by the s-t walk, as a dict that maps an edge to its multiplicity
		***
		must include s and t in the dict
		***
	path_constraints_next_path: list of subpath constraints to be satisfied, 
		each as a list of alternating nodes and edges
		***
		Note: the traversal may not satisfy all subpath constraints
		in case not all subpath constraints are satisfied, return the eulerian traversal satisfying the
		maximum number of subpath constraints
		***
	path_constraints_support: num long reads supporting each subpath constraint
	"""
	lseg = len(g.sequence_edges)
	endnode_list = [node for node in g.endnodes.keys()]

	eulerian_path = [] # A path is edge - node list starting and ending with edges
				# Since Eulerian, there could be subcycles in the middle of a path
	eulerian_path_ = [] # Path in AA cycle format
	best_path = [] # Path in AA cycle format
	valid = 0
	num_trials = 0
	l = len(path_constraints_next_path)
	unsatisfied_path_metric = [range(l), 100 * l, 100 * max(path_constraints_support + [0])]
	while valid <= 0 and num_trials < 1000:
		valid = 1
		num_trials += 1
		eulerian_path = [] 
		eulerian_path_ = []
		edges_cur = edges_next_path.copy()
		src_edge = ()
		last_seq_edge = lseg
		last_edge_dir = '+'
		for edge in edges_cur.keys(): #Start with the edge with smallest index
			if (edge[0] == 's' or edge[0] == 't'):
				src_edge = edge
				node = (g.source_edges[edge[1]][3], g.source_edges[edge[1]][4], g.source_edges[edge[1]][5])
				if len(eulerian_path) == 0:
					last_edge_dir = global_names.neg_plus_minus[node[2]]
					eulerian_path.append(('$', -1))
					eulerian_path.append(node)
					last_seq_edge = g.nodes[node][0][0]
				elif g.nodes[node][0][0] < last_seq_edge:
					last_edge_dir = global_names.neg_plus_minus[node[2]]
					eulerian_path[-1] = node
					last_seq_edge = g.nodes[node][0][0]
			elif (edge[0] == 'ns' or edge[0] == 'nt'):
				src_edge = edge
				node = endnode_list[edge[1]]
				if len(eulerian_path) == 0:
					last_edge_dir = global_names.neg_plus_minus[node[2]]
					eulerian_path.append(('$', -1))
					eulerian_path.append(node)
					last_seq_edge = g.nodes[node][0][0]
				elif g.nodes[node][0][0] < last_seq_edge:
					last_edge_dir = global_names.neg_plus_minus[node[2]]
					eulerian_path[-1] = node
					last_seq_edge = g.nodes[node][0][0]					
		del edges_cur[src_edge]
		eulerian_path.append(('s', last_seq_edge))
		if last_edge_dir == '+':
			eulerian_path_.append(str(last_seq_edge + 1) + '+')
		else:
			eulerian_path_.append(str(last_seq_edge + 1) + '-')
		edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
		if edges_cur[('e', last_seq_edge)] == 0:
			del edges_cur[('e', last_seq_edge)]
		while len(edges_cur) > 0:
			seq_edge = g.sequence_edges[last_seq_edge]
			node = (seq_edge[0], seq_edge[2], '+')
			if last_edge_dir == '-':
				node = (seq_edge[0], seq_edge[1], '-')
			eulerian_path.append(node)
			if len(edges_cur) == 1 and (list(edges_cur.keys())[0][0] == 's' or list(edges_cur.keys())[0][0] == 'ns' or \
				list(edges_cur.keys())[0][0] == 't' or list(edges_cur.keys())[0][0] == 'nt'):
				eulerian_path.append(('$', -1))
				break
			next_bp_edges = [] # Since cycle, only consider discordant edges and concordant edges
			for ci in g.nodes[node][1]:
				next_bp_edges.append(('c', ci))
			for di in g.nodes[node][2]:
				next_bp_edges.append(('d', di))
			del_list = [i for i in range(len(next_bp_edges)) if next_bp_edges[i] not in edges_cur]
			for i in del_list[::-1]:
				del next_bp_edges[i]
			if len(next_bp_edges) == 0:
				valid = 0
				break
			if len(next_bp_edges) == 1: # No branching on the path
				eulerian_path.append(next_bp_edges[0])
				edges_cur[next_bp_edges[0]] = int(edges_cur[next_bp_edges[0]]) - 1
				if edges_cur[next_bp_edges[0]] == 0:
					del edges_cur[next_bp_edges[0]]
				bp_edge = []
				if next_bp_edges[0][0] == 'c':
					bp_edge = g.concordant_edges[next_bp_edges[0][1]][:6]
				else:
					bp_edge = g.discordant_edges[next_bp_edges[0][1]][:6]
				node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
				if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
					node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
				eulerian_path.append(node_)
				last_seq_edge = g.nodes[node_][0][0]
				eulerian_path.append(('s', last_seq_edge))
				if node_[2] == '-':
					last_edge_dir = '+'
					eulerian_path_.append(str(last_seq_edge + 1) + '+')
				else:
					last_edge_dir = '-'
					eulerian_path_.append(str(last_seq_edge + 1) + '-')
				edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
				if edges_cur[('e', last_seq_edge)] == 0:
					del edges_cur[('e', last_seq_edge)]
			else:
				r = random.randint(0, len(next_bp_edges) - 1)
				eulerian_path.append(next_bp_edges[r])
				edges_cur[next_bp_edges[r]] = int(edges_cur[next_bp_edges[r]]) - 1
				if edges_cur[next_bp_edges[r]] == 0:
					del edges_cur[next_bp_edges[r]]
				bp_edge = []
				if next_bp_edges[r][0] == 'c':
					bp_edge = g.concordant_edges[next_bp_edges[r][1]][:6]
				else:
					bp_edge = g.discordant_edges[next_bp_edges[r][1]][:6]
				node_ = (bp_edge[0], bp_edge[1], bp_edge[2])
				if node == (bp_edge[0], bp_edge[1], bp_edge[2]):
					node_ = (bp_edge[3], bp_edge[4], bp_edge[5])
				eulerian_path.append(node_)
				last_seq_edge = g.nodes[node_][0][0]
				eulerian_path.append(('s', last_seq_edge))
				if node_[2] == '-':
					last_edge_dir = '+'
					eulerian_path_.append(str(last_seq_edge + 1) + '+')
				else:
					last_edge_dir = '-'
					eulerian_path_.append(str(last_seq_edge + 1) + '-')
				edges_cur[('e', last_seq_edge)] = int(edges_cur[('e', last_seq_edge)]) - 1
				if edges_cur[('e', last_seq_edge)] == 0:
					del edges_cur[('e', last_seq_edge)]
		if valid == 1 and len(best_path) == 0:
			best_path = eulerian_path_
		path_metric = [[], 0, 0]
		# check if the remaining path constraints are satisfied
		for pathi in range(len(path_constraints_next_path)):
			path_ = path_constraints_next_path[pathi]
			s = 0
			for ei in range(2, len(eulerian_path) - 1 - len(path_)):
				if eulerian_path[ei: ei + len(path_)] == path_[:] or eulerian_path[ei: ei + len(path_)] == path_[::-1]:
					s = 1
					break
			if s == 0 and valid == 1:
				path_metric[0].append(pathi)
				path_metric[1] += len(path_)
				path_metric[2] += path_constraints_support[pathi]
		if valid == 1 and len(path_metric[0]) > 0:
			valid = -1
		if valid != 0 and (len(path_metric[0]) < len(unsatisfied_path_metric[0])) or \
			(len(path_metric[0]) == len(unsatisfied_path_metric[0]) and path_metric[1] < unsatisfied_path_metric[1]) or \
			(len(path_metric[0]) == len(unsatisfied_path_metric[0]) and path_metric[1] == unsatisfied_path_metric[1] and \
			path_metric[2] < unsatisfied_path_metric[2]):
			unsatisfied_path_metric[0] = path_metric[0]
			unsatisfied_path_metric[1] = path_metric[1]
			unsatisfied_path_metric[1] = path_metric[1]
			best_path = eulerian_path_
	if len(unsatisfied_path_metric[0]) == 0:
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tPath satisfies all subpath constraints.")
	else:
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tThe following path constraints are not satisfied:")
		for pathi in unsatisfied_path_metric[0]:
			logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\t%s" %path_constraints_next_path[pathi])
	return best_path


