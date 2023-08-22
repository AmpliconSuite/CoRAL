import argparse


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "Convert cycle files in AA format to bed format.")
	parser.add_argument("--cycle_fn", help = "Input AA-formatted cycle file.", required = True)
	parser.add_argument("--output_fn", help = "Output file name.", required = True)
	parser.add_argument("--num_cycles", help = "If specified, only convert the first NUM_CYCLES cycles.", type = int)
	args = parser.parse_args()

	all_segs = dict()
	cycles = dict()
	with open(args.cycle_fn) as fp:
		for line in fp:
			t = line.strip().split()
			if t[0] == "Segment":
				all_segs[t[1]] = [t[2], int(t[3]), int(t[4])]
			if t[0][:5] == "Cycle":
				cycle_id = t[0].split('=')[1].split(';')[0]
				cycle_weight = float(t[0].split('=')[2].split(';')[0])
				cycle_segs = t[0].split('=')[3].split(';')[0].split(',')
				iscyclic = (cycle_segs[0] != "0+" or cycle_segs[-1] != "0-")
				cycle = []
				for seg in cycle_segs:
					segi = seg[:-1]
					segdir = seg[-1]
					if int(segi) > 0:
						if cycle == []:
							cycle.append(all_segs[segi] + [segdir])
						else:
							if cycle[-1][-1] == '+' and segdir == '+' and cycle[-1][0] == all_segs[segi][0] and cycle[-1][2] + 1 == all_segs[segi][1]:
								cycle[-1][2] = all_segs[segi][2]
							elif cycle[-1][-1] == '-' and segdir == '-' and cycle[-1][0] == all_segs[segi][0] and cycle[-1][1] - 1 == all_segs[segi][2]:
								cycle[-1][1] = all_segs[segi][1]
							else:
								cycle.append(all_segs[segi] + [segdir])
				if cycle[-1][-1] == '+' and cycle[0][-1] == '+' and cycle[-1][0] == cycle[0][0] and cycle[-1][2] + 1 == cycle[0][1]:
					cycle[0][1] = cycle[-1][1]
					del cycle[-1]
				if cycle[-1][-1] == '-' and cycle[0][-1] == '+' and cycle[-1][0] == cycle[0][0] and cycle[-1][1] - 1 == cycle[0][2]:
					cycle[0][2] = cycle[-1][2]
					del cycle[-1]
				cycles[int(cycle_id)] = [iscyclic, cycle_weight, cycle]
	
	with open(args.output_fn, 'w') as fp:
		fp.write("#chr\tstart\tend\torientation\tcycle_id\tiscyclic\tweight\n")
		num_cycles = len(cycles)
		if args.num_cycles:
			num_cycles = min(num_cycles, args.num_cycles)
		for i in range(1, num_cycles + 1):
			for seg in cycles[i][2]:
				fp.write("%s\t%d\t%d\t%s\t%d\t%s\t%f\n" %(seg[0], seg[1], seg[2], seg[3], i, cycles[i][0], cycles[i][1]))


	
