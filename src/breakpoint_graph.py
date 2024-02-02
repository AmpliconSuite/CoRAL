"""
Implements a class for BreakpointGraph. Will serve as a container for the 
sequence edges, corcordant edges, discordant edges, and breakpoints inferred
from a given sequencing file. 
"""
import time
import logging
import warnings
import numpy as np
import cvxopt
import cvxopt.modeling

import global_names


class BreakpointGraph:
	"""
	A container object for the breakpoint graphs.
	"""
	amplicon_intervals = [] # Amplicon intervals
	sequence_edges = [] # Sequence edges
	concordant_edges = [] # Concordant edges
	discordant_edges = [] # Discordant edges
	source_edges = [] # Source edges

	"""
	nodes: adjacent list - keys with format (chr, pos, orientation); 
	vals = [[sequence edges], [concordant edges], [discordant edges], [source edges]]  
	"""
	nodes = dict()
	endnodes = dict()
	max_cn = 0.0


	def __init__(self):
		self.amplicon_intervals = []
		self.sequence_edges = [] 
		self.concordant_edges = []
		self.discordant_edges = []
		self.source_edges = []
		self.nodes = dict()
		self.endnodes = dict()
		self.max_cn = 0.0


	def add_node(self, node_):
		"""
		Add a new node to the breakpoint graph.

		Args:
			node_: A breakpoint node.
		"""
		if type(node_) != tuple or len(node_) != 3:
			raise Exception("Breakpoint node must be of form (chr, pos, orientation).")
		if node_ in self.nodes:
			pass
		self.nodes[node_] = [[], [], [], []]


	def add_endnode(self, node_):
		"""
		Add a new node to the list corresponding to interval ends.

		Args:
			node_: A breakpoint node.
		"""
		if type(node_) != tuple or len(node_) != 3:
			raise Exception("Breakpoint node must be of form (chr, pos, orientation).")
		if node_ not in self.endnodes:
			self.endnodes[node_] = []
		else:
			warnings.warn("Node corresponding to interval end already exists.")

	
	def del_endnode(self, node_):
		"""
		Delete a node corresponding to interval ends.

		Args:
			node_: A breakpoint node.
		"""
		if node_ in self.endnodes:
			del self.endnodes[node_]
		else:
			warnings.warn("Node corresponding to interval end not exists.")


	def del_discordant_endnodes(self):
		"""
		Delete nodes correspond to interval ends and connect to a discordant edges.
		"""
		del_list = []
		for node in self.endnodes:
			if len(self.endnodes[node]) > 0:
				del_list.append(node)
		for node in del_list:
			del self.endnodes[node]


	def add_sequence_edge(self, chr, l, r, sr_count = -1, sr_flag = 'd', lr_count = -1, lr_nc = 0, cn = 0.0):
		"""
		Add a sequence edge to the graph.
		"""
		if (chr, l, '-') not in self.nodes or (chr, r, '+') not in self.nodes:
			raise Exception("Breakpoint node must be added first.")
		lseq = len(self.sequence_edges)
		self.nodes[(chr, l, '-')][0].append(lseq)
		self.nodes[(chr, r, '+')][0].append(lseq)
		self.sequence_edges.append([chr, l, r, sr_count, sr_flag, lr_count, lr_nc, r - l + 1, cn])


	def add_concordant_edge(self, chr1, pos1, o1, chr2, pos2, o2, sr_count = -1, sr_flag = 'd', lr_count = -1, reads = set([]), cn = 0.0):
		"""
		Add a concordant edge to the graph.
		"""
		if chr1 != chr2 or pos2 != pos1 + 1 or o1 != '+' or o2 != '-':
			raise Exception("Invalid concordant edge.")
		if (chr1, pos1, o1) not in self.nodes or (chr2, pos2, o2) not in self.nodes:
			raise Exception("Breakpoint node must be added first.")
		lc = len(self.concordant_edges)
		self.nodes[(chr1, pos1, o1)][1].append(lc)
		self.nodes[(chr2, pos2, o2)][1].append(lc)
		self.concordant_edges.append([chr1, pos1, o1, chr2, pos2, o2, sr_count, sr_flag, lr_count, reads, cn])


	def add_discordant_edge(self, chr1, pos1, o1, chr2, pos2, o2, sr_count = -1, sr_flag = 'd', \
				sr_cn = 0.0, lr_count = -1, reads = set([]), cn = 0.0):
		"""
		Add a discordant edge to the graph.
		"""
		if (chr1, pos1, o1) not in self.nodes or (chr2, pos2, o2) not in self.nodes:
			raise Exception("Breakpoint node must be added first.")
		ld = len(self.discordant_edges)
		self.nodes[(chr1, pos1, o1)][2].append(ld)
		self.nodes[(chr2, pos2, o2)][2].append(ld)
		if (chr1, pos1, o1) in self.endnodes:
			self.endnodes[(chr1, pos1, o1)].append(ld)
		if (chr2, pos2, o2) in self.endnodes:
			self.endnodes[(chr2, pos2, o2)].append(ld)
		self.discordant_edges.append([chr1, pos1, o1, chr2, pos2, o2, sr_count, sr_flag, sr_cn, lr_count, reads, cn])

	
	def del_discordant_edges(self, del_list, bpi_map):
		"""
		Delete a list discordant edges from the graph.
		"""
		sorted_del_list = sorted(del_list, reverse = True)
		for bpi in sorted_del_list:
			del self.discordant_edges[bpi]
		for node in self.endnodes.keys():
			for i in range(len(self.endnodes[node])):
				if self.endnodes[node][i] in sorted_del_list:
					del self.endnodes[node][i]
				else:
					self.endnodes[node][i] = bpi_map[self.endnodes[node][i]]
		for node in self.nodes.keys():
			for i in range(len(self.nodes[node][2])):
				if self.nodes[node][2][i] in sorted_del_list:
					del self.nodes[node][2][i]
				else:
					self.nodes[node][2][i] = bpi_map[self.nodes[node][2][i]]


	def add_source_edge(self, chr1, pos1, o1, sr_count = 0, sr_flag = 'd', sr_cn = 0.0, lr_cn = 0.0, cn = 0.0):
		"""
		Adds a source edge to the graph.
		"""
		if (chr1, pos1, o1) not in self.nodes:
			raise Exception("Breakpoint node must be added first.")
		self.nodes[(chr1, pos1, o1)][3].append(len(self.source_edges))
		self.source_edges.append(['source', -1, '-', chr1, pos1, o1, sr_count, sr_flag, sr_cn, lr_cn, cn])


	def del_source_edges(self, del_list, srci_map):
		"""
		Delete a list source edges from the graph.
		"""
		sorted_del_list = sorted(del_list, reverse = True)
		for srci in sorted_del_list:
			del self.source_edges[srci]
		for node in self.nodes.keys():
			for i in range(len(self.nodes[node][3])):
				if self.nodes[node][3][i] in sorted_del_list:
					del self.nodes[node][3][i]
				else:
					self.nodes[node][3][i] = srci_map[self.nodes[node][3][i]]


	def del_redundant_sequence_edges(self):
		"""
		Delete redundant sequence edges after merging.
		"""
		if len(self.discordant_edges) == 0:
			return
		del_list = []
		for seqi in range(len(self.sequence_edges)):
			sseg = self.sequence_edges[seqi]
			node1 = (sseg[0], sseg[1], '-')
			node2 = (sseg[0], sseg[2], '+')
			s1 = len(self.nodes[node1][1]) + len(self.nodes[node1][2]) + len(self.nodes[node1][3])
			s2 = len(self.nodes[node2][1]) + len(self.nodes[node2][2]) + len(self.nodes[node2][3])	
			if s1 + s2 == 0:
				del_list.append(seqi)
		for seqi in del_list[::-1]:
			ai = self.sequence_edges[seqi][:3]
			if ai in self.amplicon_intervals:
				del self.amplicon_intervals[self.amplicon_intervals.index(ai)]
			node1 = (ai[0], ai[1], '-')
			node2 = (ai[0], ai[2], '+')
			del self.self.sequence_edges[seqi]
			del self.nodes[node1]
			del self.nodes[node2]
			self.del_endnode(node1)
			self.del_endnode(node2)
		for seqi in range(len(self.self.sequence_edges)):
			sseg = self.sequence_edges[seqi]
			node1 = (sseg[0], sseg[1], '-')
			node2 = (sseg[0], sseg[2], '+')
			self.nodes[node1][0][0] = seqi
			self.nodes[node2][0][0] = seqi


	def merge_edges(self):
		"""
		Merge sequence edges connected only by concordant edges;
		Delete the nodes and concordant edges accordingly.
		"""
		c_del_list = [] # The list of concordant edges to be deleted
		seq_del_list = [] # The list of sequence edges to be deleted
		for ci in range(len(self.concordant_edges)):
			ce = self.concordant_edges[ci]
			node1 = (ce[0], ce[1], ce[2])
			node2 = (ce[3], ce[4], ce[5])
			if len(self.nodes[node1][2]) == 0 and len(self.nodes[node2][2]) == 0 and \
				len(self.nodes[node1][3]) == 0 and len(self.nodes[node2][3]) == 0:
				seqi1 = self.nodes[node1][0][0]
				seqi2 = self.nodes[node2][0][0]
				seq_del_list.append(seqi1)
				del self.nodes[node1]
				del self.nodes[node2]
				c_del_list.append(ci)
		seq_del_list = sorted(seq_del_list)
		si = 0
		li = 0
		for i in range(1, len(seq_del_list)):
			if seq_del_list[i] == seq_del_list[li] + 1:
				li += 1
			else:
				seqi1 = seq_del_list[si]
				seqi2 = seq_del_list[li] + 1
				self.sequence_edges[seqi2][1] = self.sequence_edges[seqi1][1]
				self.sequence_edges[seqi2][3] = -1
				self.sequence_edges[seqi2][4] = 'f'
				self.sequence_edges[seqi2][-2] = self.sequence_edges[seqi2][2] - self.sequence_edges[seqi2][1] + 1
				si = i
				li = i
		seqi1 = seq_del_list[si]
		seqi2 = seq_del_list[li] + 1
		self.sequence_edges[seqi2][1] = self.sequence_edges[seqi1][1]
		self.sequence_edges[seqi2][3] = -1
		self.sequence_edges[seqi2][4] = 'f'
		self.sequence_edges[seqi2][-2] = self.sequence_edges[seqi2][2] - self.sequence_edges[seqi2][1] + 1
		for seqi in seq_del_list[::-1]:
			del self.sequence_edges[seqi]
		for ci in sorted(c_del_list, reverse = True):
			del self.concordant_edges[ci]
		for seqi in range(len(self.sequence_edges)):
			sseg = self.sequence_edges[seqi]
			node1 = (sseg[0], sseg[1], '-')
			node2 = (sseg[0], sseg[2], '+')
			self.nodes[node1][0][0] = seqi
			self.nodes[node2][0][0] = seqi
		for ci in range(len(self.concordant_edges)):
			ce = self.concordant_edges[ci]
			node1 = (ce[0], ce[1], ce[2])
			node2 = (ce[3], ce[4], ce[5])
			self.nodes[node1][1][0] = ci
			self.nodes[node2][1][0] = ci


	def sort_edges(self):
		"""
		Sort sequence and concordant edges according to chromosome and position
		Reset adjacent list
		"""
		self.sequence_edges = sorted(self.sequence_edges, key = lambda sseg: (global_names.chr_idx[sseg[0]], sseg[1]))
		self.concordant_edges = sorted(self.concordant_edges, key = lambda ce: (global_names.chr_idx[ce[0]], ce[1]))
		
		for seqi in range(len(self.sequence_edges)):
			sseg = self.sequence_edges[seqi]
			self.nodes[(sseg[0], sseg[1], '-')][0] = [seqi]
			self.nodes[(sseg[0], sseg[2], '+')][0] = [seqi]
		for ci in range(len(self.concordant_edges)):
			ce = self.concordant_edges[ci]
			self.nodes[(ce[0], ce[1], ce[2])][1] = [ci]
			self.nodes[(ce[3], ce[4], ce[5])][1] = [ci]


	def compute_cn_sr_lr(self, normal_cov_sr, sr_length, normal_cov_lr, downsample_factor, min_sr_alignment_length = 30):
		"""
		Estimate CN for each edge, with both short and long reads
		"""
		lseq = len(self.sequence_edges)
		lc = len(self.concordant_edges)
		ld = len(self.discordant_edges)
		lsrc = len(self.source_edges)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Adjacent list for estimating CN:")
		for node in self.nodes.keys():
			logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Node %s; adjacent list = %s." %(str(node), self.nodes[node]))
		nvariables = lseg + lc + ld + lsrc
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "There are %d variables for cvxopt." %(nvariables))
		self.del_discordant_endnodes()
		nconstraints = len([node for node in self.nodes.keys() if node not in self.endnodes])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "There are %d constraints for cvxopt." %(nconstraints))

		wcn = [(normal_cov_sr * se[7] / sr_length + 0.5 * normal_cov_lr * se[7]) for se in self.sequence_edges]
		wcn += [normal_cov_sr * (sr_length - 1.0) / sr_length + normal_cov_lr for eci in range(lc)]
		wcn += [normal_cov_sr * (sr_length - 2 * min_sr_alignment_length) / \
				sr_length + normal_cov_lr for edi in range(ld)]
		wcn += [normal_cov_sr * (sr_length - 2 * min_sr_alignment_length) / sr_length for srci in range(lsrc)]
		wlncn = []
		for se in self.sequence_edges:
			if se[4] == 'd':
				wlncn.append(se[3] * downsample_factor - 0.5)
			else:
				wlncn.append(se[3] - 0.5)
		for ce in self.concordant_edges:
			if ce[7] == 'd':
				wlncn.append(ce[6] * downsample_factor + ce[8])
			else:
				wlncn.append((ce[6] + ce[8]) * 1.0)
		for de in self.discordant_edges:
			if de[7] == 'd':
				wlncn.append(de[6] * downsample_factor + de[9])
			else:
				wlncn.append((de[6] + de[9]) * 1.0)
		for srce in self.source_edges:
			if srce[7] == 'd':
				wlncn.append(srce[6] * downsample_factor if srce[6] >= 1 else 0.1)
			else:
				wlncn.append(srce[6] * 1.0 if srce[6] >= 1 else 0.1)
		wlrseg = [(0.5 * se[6] ** 2 / (normal_cov_lr * se[7])) for se in self.sequence_edges]
		wlrseg += [0.0 for ce in self.concordant_edges]
		wlrseg += [0.0 for de in self.discordant_edges]
		wlrseg += [0.0 for es in self.source_edges]
		wcn = cvxopt.matrix(wcn)
		wlncn = cvxopt.matrix(wlncn)
		wlrseg = cvxopt.matrix(wlrseg)
		
		cidx = 0
		balance_constraints = np.zeros([nconstraints, nvariables])
		for node in self.nodes.keys():
			if node not in self.endnodes:
				for seqi in self.nodes[node][0]:
					balance_constraints[cidx][seqi] = 1
				for ci in self.nodes[node][1]:
					balance_constraints[cidx][lseq + ci] = -1
				for di in self.nodes[node][2]:
					balance_constraints[cidx][lseq + lc + di] = -1
				for srci in self.nodes[node][3]:
					balance_constraints[cidx][lseq + lc + ld + srci] = -1
				cidx += 1
		balance_constraints = cvxopt.matrix(balance_constraints)

		# Convex optimization function required by cvxopt
		def F_normal(x = None, z = None):
			if x is None: 
				return 0, cvxopt.matrix(1.0, (nvariables, 1))
			if min(x) <= 0.0: 
				return None
			f = cvxopt.modeling.dot(wlrseg, x**-1) + cvxopt.modeling.dot(wcn, x) - cvxopt.modeling.dot(wlncn, cvxopt.log(x))
			Df = (wcn - cvxopt.mul(wlncn, x**-1) - cvxopt.mul(wlrseg, x**-2)).T
			if z is None: 
				return f, Df
			H = cvxopt.spdiag(z[0] * (cvxopt.mul(wlncn, x**-2) + cvxopt.mul(2.0 * wlrseg, x**-3)))
			return f, Df, H
		
		options = {'maxiters': 1000, 'show_progress': False}
		sol = ''
		if nconstraints > 0:
			sol = cvxopt.solvers.cp(F_normal, A = balance_constraints, 
						b = cvxopt.matrix([0.0 for i in range(nconstraints)]), 
						kktsolver = 'ldl', options = options)
			if sol['status'] == "optimal" or sol['status'] == "unknown":
				if sol['status'] == "optimal":
					logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Found optimal solution.")
				else:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Reached maximum num iterations.")
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tprimal objective = %f" %(sol['primal objective']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tdual objective = %f" %(sol['dual objective']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tgap = %f" %(sol['gap']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\trelative gap = %f" %(sol['relative gap']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tprimal infeasibility = %f" %(sol['primal infeasibility']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tdual infeasibility = %f" %(sol['dual infeasibility']))
				for seqi in range(lseq):
					self.sequence_edges[seqi][-1] = sol['x'][seqi] * 2
					if sol['x'][seqi] * 2 > self.max_cn:
						self.max_cn = sol['x'][seqi] * 2
				for ci in range(lc):
					self.concordant_edges[ci][-1] = sol['x'][lseq + ci] * 2
					if sol['x'][lseq + ci] * 2 > self.max_cn:
						self.max_cn = sol['x'][lseq + ci] * 2
				for di in range(ld):
					self.discordant_edges[di][-1] = sol['x'][lseq + lc + di] * 2
					if sol['x'][lseq + lc + di] * 2 > self.max_cn:
						self.max_cn = sol['x'][lseq + lc + di] * 2
				for srci in range(len(self.source_edges)):
					self.source_edges[srci][-1] = sol['x'][lseq + lc + ld + srci] * 2
					if sol['x'][lseq + lc + ld + srci] * 2 > self.max_cn:
						self.max_cn = sol['x'][lseq + lc + ld + srci] * 2
		else:
			assert lc == 0 and ld == 0 and lsrc == 0
			logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Skipped convex optimization.")
			for seqi in range(lseq):
				se = self.sequence_edges[seqi]
				cn_seqi = 0.0
				if se[4] == 'd':
					cn_seqi = (sr_length * se[3]) / (10.0 * se[7])
				else:
					cn_seqi = (sr_length * se[3]) / (normal_cov_sr * se[7])
				cn_seqi += se[6] / (normal_cov_lr * se[7]) 
				self.sequence_edges[seqi][-1] = cn_seqi
				if cn_seqi > self.max_cn:
					self.max_cn = cn_seqi
		self.max_cn += 1.0


	def compute_cn_lr(self, normal_cov_lr):
		"""
		Estimate CN for each edge, with only long reads
		"""
		lseq = len(self.sequence_edges)
		lc = len(self.concordant_edges)
		ld = len(self.discordant_edges)
		lsrc = len(self.source_edges)
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Adjacent list for estimating CN:")
		for node in self.nodes.keys():
			logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Node %s; adjacent list = %s." %(str(node), self.nodes[node]))
		nvariables = lseq + lc + ld + lsrc
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "There are %d variables for cvxopt." %(nvariables))
		nconstraints = len([node for node in self.nodes.keys() if node not in self.endnodes])
		logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "There are %d constraints for cvxopt." %(nconstraints))

		wcn = []
		wlncn = []
		wlrseg = []
		wcn = [0.5 * normal_cov_lr * se[7] for se in self.sequence_edges]
		wcn += [normal_cov_lr for eci in range(lc)]
		wcn += [normal_cov_lr for edi in range(ld)]
		wcn += [0.5 * normal_cov_lr for srci in range(lsrc)]
		wlncn = [-0.5 for seg in self.sequence_edges]
		wlncn += [ce[8] * 1.0 for ce in self.concordant_edges]
		wlncn += [de[9] * 1.0 for de in self.discordant_edges]
		wlncn += [-0.5 for srci in range(lsrc)]
		wlrseg = [(0.5 * se[6] ** 2 / (normal_cov_lr * se[7])) for se in self.sequence_edges]
		wlrseg += [0.0 for eci in range(lc)]
		wlrseg += [0.0 for edi in range(ld)]
		wlrseg += [(0.5 * self.source_edges[srci][-1] ** 2 / normal_cov_lr) for srci in range(lsrc)]
		wcn = cvxopt.matrix(wcn)
		wlncn = cvxopt.matrix(wlncn)
		wlrseg = cvxopt.matrix(wlrseg)
		
		cidx = 0
		balance_constraints = np.zeros([nconstraints, nvariables])
		for node in self.nodes.keys():
			if node not in self.endnodes:
				for seqi in self.nodes[node][0]:
					balance_constraints[cidx][seqi] = 1
				for eci in self.nodes[node][1]:
					balance_constraints[cidx][lseq + eci] = -1
				for edi in self.nodes[node][2]:
					balance_constraints[cidx][lseq + lc + edi] = -1
				for srci in self.nodes[node][3]:
					balance_constraints[cidx][lseq + lc + ld + srci] = -1
				cidx += 1
		balance_constraints = cvxopt.matrix(balance_constraints)

		# Convex optimization function required by cvxopt
		def F_normal(x = None, z = None):
			if x is None: 
				return 0, cvxopt.matrix(1.0, (nvariables, 1))
			if min(x) <= 0.0: 
				return None
			f = cvxopt.modeling.dot(wlrseg, x**-1) + cvxopt.modeling.dot(wcn, x) - cvxopt.modeling.dot(wlncn, cvxopt.log(x))
			Df = (wcn - cvxopt.mul(wlncn, x**-1) - cvxopt.mul(wlrseg, x**-2)).T
			if z is None: 
				return f, Df
			H = cvxopt.spdiag(z[0] * (cvxopt.mul(wlncn, x**-2) + cvxopt.mul(2.0 * wlrseg, x**-3)))
			return f, Df, H
				
		options = {'maxiters': 1000, 'show_progress': False}
		sol = ''
		if nconstraints > 0:
			sol = cvxopt.solvers.cp(F_normal, A = balance_constraints, 
						b = cvxopt.matrix([0.0 for i in range(nconstraints)]), 
						kktsolver = 'ldl', options = options)
			if sol['status'] == "optimal" or sol['status'] == "unknown":
				if sol['status'] == "optimal":
					logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Found optimal solution.")
				else:
					logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Reached maximum num iterations.")
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tprimal objective = %f" %(sol['primal objective']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tdual objective = %f" %(sol['dual objective']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tgap = %f" %(sol['gap']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\trelative gap = %f" %(sol['relative gap']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tprimal infeasibility = %f" %(sol['primal infeasibility']))
				logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "\tdual infeasibility = %f" %(sol['dual infeasibility']))
				for seqi in range(lseq):
					self.sequence_edges[seqi][-1] = sol['x'][seqi] * 2
					if sol['x'][seqi] * 2 > self.max_cn:
						self.max_cn = sol['x'][seqi] * 2
				for ci in range(lc):
					self.concordant_edges[ci][-1] = sol['x'][lseq + ci] * 2
					if sol['x'][lseq + ci] * 2 > self.max_cn:
						self.max_cn = sol['x'][lseq + ci] * 2
				for di in range(ld):
					de = self.discordant_edges[di]
					if de[0] == de[3] and de[1] == de[4] and de[2] == de[5]:
						self.discordant_edges[di][-1] = sol['x'][lseq + lc + di]
						if sol['x'][lseq + lc + di] > self.max_cn:
							self.max_cn = sol['x'][lseq + lc + di]
					else:
						self.discordant_edges[di][-1] = sol['x'][lseq + lc + di] * 2
						if sol['x'][lseq + lc + di] * 2 > self.max_cn:
							self.max_cn = sol['x'][lseq + lc + di] * 2
				for srci in range(lsrc):
					self.source_edges[srci][-1] = sol['x'][lseq + lc + ld + srci] * 2
					if sol['x'][lseq + lc + ld + srci] * 2 > self.max_cn:
						self.max_cn = sol['x'][lseq + lc + ld + srci] * 2
		else:
			assert lc == 0 and ld == 0 and lsrc == 0
			logging.debug("#TIME " + '%.4f\t' %(time.time() - global_names.TSTART) + "Skipped convex optimization.")
			for seqi in range(lseq):
				se = self.sequence_edges[seqi]
				cn_seqi = se[6] * 2.0 / (normal_cov_lr * se[7])
				self.sequence_edges[seqi][-1] = cn_seqi
				if cn_seqi > self.max_cn:
					self.max_cn = cn_seqi
		self.max_cn += 1.0


	def infer_max_seq_repeat(self, gain = 5.0, size_cutoff = 10000, repeat = 2):
		"""
		Estimate maximum allowed multiplicities in cycles/paths for any sequence edge

		gain: float, only consider sequence edges with predicted CN >= gain
		size_cutoff: integer, only consider sequence edges of size CN >= size_cutoff
		repeat: integer, default maximum multiplicity
			default value is 2, indicating both orientations

		Return: integer, estimated maximum multiplicity on sequence edges
		"""
		max_seq_repeat = repeat
		seq_cn_list = [se[-1] for se in self.sequence_edges if se[7] >= size_cutoff and se[-1] >= gain]
		seq_len_list = [se[7] for se in self.sequence_edges if se[7] >= size_cutoff and se[-1] >= gain]
		if len(seq_cn_list) > 0:
			max_cn = max(seq_cn_list)
			avg_cn = np.average(seq_cn_list, weights = seq_len_list)
			max_seq_repeat = int(round(max_cn / avg_cn)) + 1
		return max_seq_repeat


	def nextminus(self, chr, pos, min_bp_match_cutoff_ = 100):
		"""
		Helper function to read_graph
		Return the distance to the next position towards 3' which has incoming breakpoint edges on chr 
		"""
		cr = -1
		pos_ = pos
		while (chr, pos_, '-') in self.nodes.keys():
			if pos_ != pos and len(self.nodes[(chr, pos_, '-')][2]) > 0:
				break
			if cr >= min_bp_match_cutoff_:
				break
			seglen = self.sequence_edges[self.nodes[(chr, pos_, '-')][0][0]][7]
			cr = max(cr, 0) + seglen
			pos_ = pos_ + seglen
		return cr


	def lastminus(self, chr, pos, min_bp_match_cutoff_ = 100):
		"""
		Helper function to read_graph
		Return the distance to the next position towards 5' which has incoming breakpoint edges on chr 
		"""
		cl = -1
		pos_ = pos
		while (chr, pos_ - 1, '+') in self.nodes.keys():
			if pos_ != pos and len(self.nodes[(chr, pos_, '-')][2]) > 0:
				break
			if cl >= min_bp_match_cutoff_:
				break
			seglen = self.sequence_edges[self.nodes[(chr, pos_ - 1, '+')][0][0]][7]
			cl = max(cl, 0) + seglen
			pos_ = pos_ - seglen
		return cl


	def nextplus(self, chr, pos, min_bp_match_cutoff_ = 100):
		"""
		Helper function to read_graph
		Return the next position towards 3' which has outgoing breakpoint edges on chr 
		"""
		cr = -1
		pos_ = pos
		while (chr, pos_ + 1, '-') in self.nodes.keys():
			if pos_ != pos and len(self.nodes[(chr, pos_, '+')][2]) > 0:
				break
			if cr >= min_bp_match_cutoff_:
				break
			seglen = self.sequence_edges[self.nodes[(chr, pos_ + 1, '-')][0][0]][7]
			cr = max(cr, 0) + seglen
			pos_ = pos_ + seglen
		return cr


	def lastplus(self, chr, pos, min_bp_match_cutoff_ = 100):
		"""
		Helper function to read_graph
		Return the next position towards 5' which has outgoing breakpoint edges on chr 
		"""
		cl = -1
		pos_ = pos
		while (chr, pos_, '+') in self.nodes.keys():
			if pos_ != pos and len(self.nodes[(chr, pos_, '+')][2]) > 0:
				break
			if cl >= min_bp_match_cutoff_:
				break
			seglen = self.sequence_edges[self.nodes[(chr, pos_, '+')][0][0]][7]
			cl = max(cl, 0) + seglen
			pos_ = pos_ - seglen
		return cl


def output_breakpoint_graph_sr_lr(g, ogfile, downsample_factor):
	"""
	Write a breakpoint graph to file in AA graph format with short read and long read information
	"""
	with open(ogfile, 'w') as fp:
		fp.write("SequenceEdge: StartPosition, EndPosition, PredictedCN, NumberOfReadPairs, NumberOfLongReads, Size\n")
		for se in g.sequence_edges:
			if se[4] == 'd':
				fp.write("sequence\t%s:%s-\t%s:%s+\t%f\t%d\t%d\t%d\n" %(se[0], se[1], se[0], se[2], se[-1],
					int(math.round(se[3] * downsample_factor)), se[5], se[7]))
			else:
				fp.write("sequence\t%s:%s-\t%s:%s+\t%f\t%d\t%d\t%d\n" 
					%(se[0], se[1], se[0], se[2], se[-1], se[3], se[5], se[7]))
		fp.write("BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfReadPairs, NumberOfLongReads\n")
		for srce in g.source_edges:
			if srce[7] == 'd':
				fp.write("source\t%s:%s%s->%s:%s%s\t%f\t-1\t%d\n" %(srce[0], srce[1], srce[2], srce[3],
					srce[4], srce[5], srce[-1], int(math.round(srce[6] * downsample_factor))))
			else:
				fp.write("source\t%s:%s%s->%s:%s%s\t%f\t-1\t%d\n" 
					%(srce[0], srce[1], srce[2], srce[3], srce[4], srce[5], srce[-1], srce[6]))
		for ce in self.concordant_edges:
			if ce[7] == 'd':
				fp.write("concordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\n" %(ce[0], ce[1], ce[2], ce[3], 
					ce[4], ce[5], ce[-1], int(math.round(ce[6] * downsample_factor)), ce[8]))
			else:
				fp.write("concordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\n" %(ce[0], ce[1], ce[2], ce[3], 
					ce[4], ce[5], ce[-1], ce[6], ce[8]))
		for de in self.discordant_edges:
			if de[7] == 'd':
				fp.write("discordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\n" %(de[0], de[1], de[2], de[3], 
					de[4], de[5], de[-1], int(math.round(de[6] * downsample_factor)), de[9]))
			else:
				fp.write("discordant\t%s:%s%s->%s:%s%s\t%f\t%d\t%d\n" %(de[0], de[1], de[2], de[3], 
					de[4], de[5], de[-1], de[6], de[9]))


def output_breakpoint_graph_lr(g, ogfile):
	"""
	Write a breakpoint graph to file in AA graph format with only long read information
	"""
	with open(ogfile, 'w') as fp:
		fp.write("SequenceEdge: StartPosition, EndPosition, PredictedCN, AverageCoverage, Size, NumberOfLongReads\n")
		for se in g.sequence_edges:
			fp.write("sequence\t%s:%s-\t%s:%s+\t%f\t%f\t%d\t%d\n" 
					%(se[0], se[1], se[0], se[2], se[-1], se[6] * 1.0 / se[7], se[7], se[5]))
		fp.write("BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfLongReads\n")
		for srce in g.source_edges:
			fp.write("source\t%s:%s%s->%s:%s%s\t%f\t-1\n" %(srce[0], srce[1], srce[2], srce[3], srce[4], srce[5], srce[-1]))
		for ce in g.concordant_edges:
			fp.write("concordant\t%s:%s%s->%s:%s%s\t%f\t%d\n" %(ce[0], ce[1], ce[2], ce[3], 
							ce[4], ce[5], ce[-1], ce[8]))
		for de in g.discordant_edges:
			fp.write("discordant\t%s:%s%s->%s:%s%s\t%f\t%d\n" %(de[0], de[1], de[2], de[3], 
				de[4], de[5], de[-1], de[9]))


def output_breakpoint_info_sr_lr(g, obpfile, downsample_factor, new_bp_stats):
	"""
	Write the list of breakpoints to file
	"""
	with open(obpfile, 'w') as fp:
		fp.write("chr1\tpos1\tchr2\tpos2\torientation\tsr_support\tlr_support\tlr_info=[avg1, avg2, std1, std2, mapq1, mapq2]\n")
		for di in range(len(g.discordant_edges)):
			de = g.discordant_edges[di]
			if di in bp_stats:
				fp.write("%s\t%s\t%s\t%s\t%s%s\t-1\t%d\t%s\n" 
				%(de[3], de[4], de[0], de[1], de[5], de[2], de[9], new_bp_stats[di]))
			else:
				if de[7] == 'd':
					fp.write("%s\t%s\t%s\t%s\t%s%s\t%d\t%d\tN/A\n" %(de[3], de[4], de[0], de[1], de[5], de[2], 
						int(math.round(de[6] * downsample_factor)), de[9]))
				else:
					fp.write("%s\t%s\t%s\t%s\t%s%s\t%d\t%d\tN/A\n" %(de[3], de[4], de[0], de[1], de[5], de[2], 
						de[6], de[9]))


def output_breakpoint_info_lr(g, obpfile, bp_stats):
	"""
	Write the list of breakpoints to file
	"""
	with open(obpfile, 'w') as fp:
		fp.write("chr1\tpos1\tchr2\tpos2\torientation\tlr_support\tlr_info=[avg1, avg2, std1, std2, mapq1, mapq2]\n")
		for di in range(len(g.discordant_edges)):
			de = g.discordant_edges[di]
			fp.write("%s\t%s\t%s\t%s\t%s%s\t%d\t%s\n" 
				%(de[3], de[4], de[0], de[1], de[5], de[2], de[9], bp_stats[di]))


