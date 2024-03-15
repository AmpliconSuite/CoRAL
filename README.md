#  CoRAL - Complete Reconstruction of Amplifications with Long reads
## Reference
https://www.biorxiv.org/content/10.1101/2024.02.15.580594v1

## Installation
CoRAL can be installed and run on most modern Unix-like operating systems (e.g. Ubuntu 18.04+, CentOS 7+, macOS). 

CoRAL requires python>=3.7, and the following packages.
* pysam(>=0.1.7) https://pysam.readthedocs.io/en/stable/ for reading mapped sequences in ```*.BAM``` format
* cvxopt https://cvxopt.org/ for estimating CN in breakpoint graph.
* Gurobi (>=9.1, for Python) https://www.gurobi.com/documentation/current/refman/py_python_api_overview.html for solving the quadratic (constrained) program to extract cycles/paths from the breakpoint graph.
* CNVkit https://cnvkit.readthedocs.io/ for producing the copy number segments, as well as seed amplification intervals, in amplified interval search.

## Command line arguments to run CoRAL
CoRAL includes a main script ```infer_breakpoint_graph.py```, which implements its main functionality of constructing breakpoint graphs and extracting cycles(/paths) from the breakpoint graph. In addition, it also includes the following 4 separate modules:
* ```cnv_seed.py``` - for seed amplification interval identification from CNV calls;
* ```plot_amplicons.py``` - for outputing stylistic visualizations of the breakpoint graphs and cycles;
* ```cycle2bed.py``` - for converting the cycles output in [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) format to ```*.bed``` format;
* ```hsr.py``` - for identifying potential integration sites on reference genome in ecDNA/HSR positive samples.

### 1. Running ```infer_breakpoint_graph.py```
Usage: 
```python3 infer_breakpoint_graph.py <Required arguments> <Optional arguments>```

**1.1 Required arguments:**
* ```--lr_bam <FILE>``` - Coordinate sorted ```*.BAM``` file, with ```*.bai``` index (mapped to the provided reference genome) in the same directory.
* ```--seed <FILE>``` - ```*.bed``` file with a putative list of seed amplification intervals. The seed amplification intervals can be obtained through [running ```cnv_seed.py```](#running-```cnv_seed.py```), or provided manually.
* ```--output_prefix <STRING>``` - Prefix of the output ```graph.txt``` and ```cycles.txt``` files. Note that if these files are desired to be written to a different directory, then paths should also be included. 
* ```--cnseg <FILE>``` - CNV segmentation (```*.cns```) file given by CNVkit.

**1.2 Optional arguments:**
* ```--min_bp_support <FLOAT>``` - Filter out breakpoints with less than (min_bp_support * normal coverage) long read support in breakpoint graph construction.
* ```--output_all_path_constraints``` - If specified, output all path constraints given by long reads in ```*.cycles``` file (see "Expected output" below).
* ```--cycle_decomp_alpha <FLOAT between [0, 1]>``` - Parameter used to balance CN weight and path constraints in the objective function of greedy cycle extraction. Default value is 0.01, higher values favor the satisfaction of more path constraints.
* ```--cycle_decomp_time_limit <INT>``` - Maximum running time (in seconds) reserved for solving the quadratic program with Gurobi (integer program solver). The solver would return the best solution(s) it currently found, regardless of the optimality status, when reaching this time limit. Default value is 7200 (i.e., 2 hours).
* ```--cycle_decomp_threads <INT>``` - Number of threads reserved for for solving the quadratic program with Gurobi (integer program solver). If not specified (and by default), the solver would attempt to use up all available cores in the working machine. 
* ```--postprocess_greedy_sol``` - If specified, automatically postprocess the cycles/paths returned in greedy cycle extraction, by solving the full quadratic program to minimize the number of cycles/paths starting with the greedy cycle extraction solution (as an initial solution).
*	```--log_fn <FILE>``` - Name of the main ```*.log``` file, which can be used to trace the status of ```infer_breakpoint_graph.py``` run(s). 

**1.3 Expected output:**

CoRAL may identify and reconstruct a few distinct focal amplifications in the input ```*.BAM``` sample, each will be organized as an *amplicon*, which includes a connected component of amplified intervals and their connections by discordant edges. CoRAL writes the following files to the path specified with ```--output_prefix```, all with the prefix given by ```--output_prefix```.

* Graph file: For each amplicon, a tab-separated text file named ```output_prefix_amplicon*_graph.txt``` describing the *sequence edges*, *concordant edges* and *discordant edges* in the graph and their predicted copy count. Note that the graph files output by CoRAL have the same format as those output by [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (and therefore the files can be used interchangeably with AmpliconArchitect). Here is an example graph file from GBM39, a cell line with *EGFR* amplified on ecDNA.
```
SequenceEdge: StartPosition, EndPosition, PredictedCN, AverageCoverage, Size, NumberOfLongReads
sequence	chr7:54659673-	chr7:54763281+	4.150534	45.907363	103609	576
sequence	chr7:54763282-	chr7:55127266+	89.340352	1052.714362	363985	40637
sequence	chr7:55127267-	chr7:55155020+	2.843655	32.729552	27754	172
sequence	chr7:55155021-	chr7:55609190+	89.340352	1013.182857	454170	49675
sequence	chr7:55609191-	chr7:55610094+	2.868261	31.027655	904	915
sequence	chr7:55610095-	chr7:56049369+	89.340352	1023.280633	439275	49106
sequence	chr7:56049370-	chr7:56149664+	4.150534	49.623899	100295	562
BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfLongReads
concordant	chr7:54763281+->chr7:54763282-	4.150534	26
concordant	chr7:55127266+->chr7:55127267-	2.843655	36
concordant	chr7:55155020+->chr7:55155021-	2.843655	32
concordant	chr7:55609190+->chr7:55609191-	2.868261	90
concordant	chr7:55610094+->chr7:55610095-	2.868261	93
concordant	chr7:56049369+->chr7:56049370-	4.150534	45
discordant	chr7:55610095-->chr7:55609190+	86.472091	869
discordant	chr7:56049369+->chr7:54763282-	85.189818	981
discordant	chr7:55155021-->chr7:55127266+	86.496697	978
```
* Cycles file: 
For each amplicon, a tab-separated text file named ```output_prefix_amplicon*_cycles.txt``` describing the list of cycles and paths returned from cycle extraction. Note that the cycles files output by CoRAL have mostly the same format as those output by [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (and therefore the files can be used interchangeably with AmpliconArchitect in most cases). Specifically a cycles file includes (i) the list of amplified intervals; (ii) the list of sequence edges; (iii) the list of cycles and paths, where an entry starts with ```0+``` and ends with ```0-``` in ```Segments``` indicates a path - these lines have the same format as AmpliconArchitect output. CoRAL's cycles files additionally include (iv) a list of longest (i.e., there are no paths that can form a sub/super-path to each other) path constraint indicated by long reads, and used in CoRAL's cycle extraction. Here is an example cycles file corresponding to the above graph file from GBM39.
```
Interval	1	chr7	54659673	56149664
List of cycle segments
Segment	1	chr7	54659673	54763281
Segment	2	chr7	54763282	55127266
Segment	3	chr7	55127267	55155020
Segment	4	chr7	55155021	55609190
Segment	5	chr7	55609191	55610094
Segment	6	chr7	55610095	56049369
Segment	7	chr7	56049370	56149664
List of longest subpath constraints
Path constraint	1	2+,3+,4+	Support<=6	Satisfied
Path constraint	2	4+,5+,6+	Support<=34	Satisfied
Cycle=1;Copy_count=82.34616279663038;Segments=2+,4+,6+;Path_constraints_satisfied=
Cycle=2;Copy_count=2.8436550275157644;Segments=0+,2+,3+,4+,5+,6+,0-;Path_constraints_satisfied=1,2
```
Note that if ```--output_all_path_constraints``` is specified, then all path constraints given by long reads will be written to in ```*.cycles``` file.
* Other outputs include the ```output_prefix_amplicon*_model.lp``` file(s) and ```output_prefix_amplicon*_model.log``` file(s) given by Gurobi (integer program solver), for each amplicon, respectively describing the quadratic (constrainted) program in a human readable format, and the standard output produced by Gurobi.
### 2. Running ```cnv_seed.py```
As the seed amplification intervals are required by the main script ```infer_breakpoint_graph.py```, it is suggested the user first run ```cnv_seed.py``` to generate seed amplification intervals.

Usage: 
```python3 cnv_seed.py <Required arguments> <Optional arguments>```

**Required arguments:**
* ```--cn_seg <FILE>```, The CNV segmentation (```*.cns```) file given by CNVkit, which is also required by ```infer_breakpoint_graph.py```.

**Optional arguments:**
* ```--out <STRING>``` - Prefix of the output ```*_CNV_SEEDS.bed``` file. Note that if these file is desired to be written to a different directory, then a path/directory should also be included. If not specified (by default), output the ```*_CNV_SEEDS.bed``` file to the current directory with the same prefix as the input ```*.cns``` file.
* ```--gain <FLOAT>``` - A minimum CN threshold (with the assumption of diploid genome) for a particular CN segment to be considered as a seed. Default is 6.0.
* ```--min_seed_size <INT>``` - Minimum size (in bp) for a CN segment to be considered as a seed. Default is 100000.
* ```--max_seg_gap <INT>``` - Maximum gap size (in bp) to merge two proximal CN segments to be considered as seed intervals. If at least two segments are merged, then they will be treated as a single candidate to be filtered with ```--min_seed_size```, and their aggregate size will be compared with the value. Default is 300000. 

### 3. Running ```cycle2bed.py```
CoRAL provides an option to convert its cycles output in AmpliconArchitect format ```*_cycles.txt``` into ```*.bed``` format (similar to [Decoil](https://github.com/madagiurgiu25/decoil-pre)), which makes it easier for downstream analysis of these cycles.

Usage: 
```python3 cycle2bed.py <Required arguments> <Optional arguments>```

**Required arguments:**
* ```--cycle_fn <FILE>``` - Input cycles file in AmpliconArchitect format.
* ```--output_fn <FILE>```  - Output cycles file in ```*.bed``` format.

**Optional arguments:** 
* ```--num_cycles <INT>``` - If specified, only convert the first NUM_CYCLES cycles.

Here is an example output of ```cycle2bed.py``` given by the above cycles file from GBM39.
```
#chr	start	end	orientation	cycle_id	iscyclic	weight
chr7	54763282	55127266	+	1	True	82.346163
chr7	55155021	55609190	+	1	True	82.346163
chr7	55610095	56049369	+	1	True	82.346163
chr7	54763282	56049369	+	2	False	2.843655
```

