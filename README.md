#  CoRAL - Complete Reconstruction of Amplifications with Long reads
## Reference
CoRAL is a tool which utilizes aligned, single-molecule long-read data (.bam) as input, and identifies candidate ecDNA structures. The original Genome Research '24 paper is available here: https://genome.cshlp.org/content/34/9/1344.

- **CoRAL only works on long-read
whole-genome sequencing data (PacBio, Oxford Nanopore, etc.) - not targeted 
sequencing!**
- **We also only support hg38-aligned data currently. Support for other genomes 
is [coming soon](https://github.com/AmpliconSuite/CoRAL/issues/33)!**

## Installation
CoRAL can be installed and run on most modern Unix-like operating systems (e.g. Ubuntu 18.04+, CentOS 7+, macOS). 

CoRAL requires python>=3.12; we recommend using venv/conda for managing Python/pip installations.

1. Clone source

    ```
    git clone https://github.com/AmpliconSuite/CoRAL
    cd CoRAL
    ```

2. Install packages using `poetry`. 
    
      ```bash
      pip install poetry
      poetry install
     ```


3. [Download a Gurobi optimizer license](https://support.gurobi.com/hc/en-us/articles/360040541251-How-do-I-obtain-a-free-academic-license) (free for academic use)
   - Place the `gurobi.lic` file you download into `$HOME/`. This path is usually `/home/username/gurobi.lic`.



4. Finish installing CNVkit dependencies (recommended)
   ```bash
   Rscript -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
   Rscript -e 'BiocManager::install("DNAcopy")'
   ```

[//]: # (* pysam&#40;>=0.1.7&#41; https://pysam.readthedocs.io/en/stable/ for reading mapped sequences in ```*.BAM``` format)

[//]: # (* cvxopt https://cvxopt.org/ for estimating CN in breakpoint graph.)

[//]: # (* Gurobi &#40;>=9.1, for Python&#41; https://www.gurobi.com/documentation/current/refman/py_python_api_overview.html for solving the quadratic &#40;constrained&#41; program to extract cycles/paths from the breakpoint graph.)

[//]: # (* CNVkit https://cnvkit.readthedocs.io/ for producing the copy number segments, as well as seed amplification intervals, in amplified interval search.)

## Getting copy number calls
Before running CoRAL, you will need genome-wide copy number (CN) calls generated from your long-read data. 

- If you have these already, simply ensure that they are in a .bed format like so:

   `chrom  start   end   CN`


- If you don't have these then you can run CNVkit (installed as a dependency) to generate them, by doing

   `./scripts/call_cnvs.sh <input.bam> ./reference/hg38full_ref_5k.cnn <output_dir>`

   This will create a file called `[input].cns`, which you can feed to CoRAL for it's `--cn_segs` argument.

## Command line arguments to run CoRAL

CoRAL and its various run-modes can by used in the following manner

`coral [mode] [mode arguments]`

The modes are as follows:
1. `seed`: Identify and filter copy number gain regions where amplifications exist
2. `reconstruct`: Perform breakpoint graph construct and cycle decomposition on the amplified seeds.
3. `plot`: Create plots of decomposed cycles and/or breakpoint graph sashimi plot.
4. `hsr`: Identify candidate locations of chromosomal homogenously staining region (HSR) integration points for ecDNA.
5. `cycle2bed`: Convert the [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (AA) style  `*_cycles.txt` file to a .bed format. The AA format is also used by CoRAL.
6. `cycle`: Run the cycle extraction algorithm on a previously generated 
breakpoint graph. NOTE: This requires the breakpoint graph to be generated with
CoRAL v2.1.0 or later, as we require `path constraints` and `amplicon intervals`
to be included in the provided `*_graph.txt` file.


## 1. ```seed```
As the seed amplification intervals are required by the main script ```reconstruct``` mode, it is suggested the user first run ```seed``` mode to generate seed amplification intervals.

Usage: 
```coral seed <Required arguments> <Optional arguments>```

**Required arguments:**
* ```--cn-seg <FILE>```, Long read segmented whole genome CN calls (.bed or CNVkit .cns file).

**Optional arguments:**
* ```--output-prefix <STRING>``` - Prefix of the output ```*_CNV_SEEDS.bed``` file.  If not specified (by default), output the ```*_CNV_SEEDS.bed``` with the same prefix as the input ```*.cns``` file.
* ```--gain <FLOAT>``` - A minimum CN threshold (with the assumption of diploid genome) for a particular CN segment to be considered as a seed. Default is 6.0.
* ```--min-seed-size <INT>``` - Minimum size (in bp) for a CN segment to be considered as a seed. Default is 100000.
* ```--max-seg-gap <INT>``` - Maximum gap size (in bp) to merge two proximal CN segments to be considered as seed intervals. If at least two segments are merged, then they will be treated as a single candidate to be filtered with ```--min-seed-size```, and their aggregate size will be compared with the value. Default is 300000. 


## 2. ```reconstruct```
Usage: 
```reconstruct <Required arguments> <Optional arguments>```

**2.1 Required arguments:**
* ```--lr-bam <FILE>``` - Coordinate sorted ```*.BAM``` file, with ```*.bai``` index (mapped to the provided reference genome) in the same directory.
* ```--cnv-seed <FILE>``` - ```*.bed``` file with a putative list of seed amplification intervals. The seed amplification intervals can be obtained through [running ```seed``` mode](#CoRAL.py-```seed```), or provided manually.
* ```--output-dir <FOLDER>``` - Directory to which the output ```graph.txt``` and ```cycles.txt``` files will be written.
* ```--cn-seg <FILE>``` - Long read segmented whole genome CN calls (.bed or CNVkit .cns file).

**2.2 Optional arguments:**
* ```--min-bp-support <FLOAT>``` - Filter out breakpoints with less than (min_bp_support * normal coverage) long read support in breakpoint graph construction. The default value is set to 1.0, meaning to filter out breakpoints supported by less than diploid coverage, but ***it is highly recommended to specify a much larger value, e.g. 10.0 to obtain a cleaner breakpoint graph and the dominating ecDNA cycle(s).***
* ```--skip-cycle-decomp``` - If specified, will stop by only outputting the breakpoint graph files ```*_graph.txt``` (see [**Expected output**](#2.3-Expected-output) below) for all amplicons and not extract cycles from the graph and output ```*_cycles.txt```.
* ```--output-all-path_constraints``` - If specified, output all path constraints given by long reads in ```*_cycles.txt``` file (see "Expected output" below).
* ```--cycle-decomp-alpha <FLOAT between [0, 1]>``` - Parameter used to balance CN weight and path constraints in the objective function of greedy cycle extraction. Default value is 0.01, higher values favor the satisfaction of more path constraints.
* ```--solver-time-limit <INT>``` - Maximum running time (in seconds) reserved for solving a single quadratic program using the chosen integer program solver (e.g. Gurobi, SCIP). The solver would return the best solution(s) it currently found, regardless of the optimality status, when reaching this time limit. Default value is 7200 (i.e., 2 hours).
* ```--solver-threads <INT>``` - Number of threads reserved for for solving the quadratic program with Gurobi (integer program solver). If not specified (and by default), the solver would attempt to use up all available cores in the working machine. 
* ```--solver <choice>``` - Solver for cycle extraction. Must be one of `[gurobi_direct, scip]`.
* ```--global-time-limit <INT>``` - Maximum running time (in seconds) reserved for the entire cycle decomposition process. Default value is 21600 (i.e., 6 hours).
* ```--postprocess-greedy-sol``` - If specified, automatically postprocess the cycles/paths returned in greedy cycle extraction, by solving the full quadratic program to minimize the number of cycles/paths starting with the greedy cycle extraction solution (as an initial solution).
*	```--log-file <FILE>``` - Name of the main ```*.log``` file, which can be used to trace the status of ```reconstruct``` run(s). 

**2.3 Expected output:**

CoRAL may identify and reconstruct a few distinct focal amplifications in the input ```*.BAM``` sample, each will be organized as an *amplicon*, which includes a connected component of amplified intervals and their connections by discordant edges. CoRAL writes the following files to the directory specified with ```--output_dir```.

* Graph file: For each amplicon, a tab-separated text file named ```output_dir/amplicon*_graph.txt``` describing the *sequence edges*, *concordant edges* and *discordant edges* in the graph and their predicted copy count. Note that the graph files outputted by CoRAL have the same format as those outputted by [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (and therefore the files can be used interchangeably with AmpliconArchitect). Here is an example graph file from GBM39, a cell line with *EGFR* amplified on ecDNA.
   * As of version 2.1.0, CoRAL additionally includes `path constraints` and 
   `amplicon intervals` in the `*_graph.txt` file. This results in the graph
   being fully self-contained and able to be passed to cycle extraction without 
   re-parsing the BAM file. For more information on how to interpret this
   metadata, visit our [wiki](https://github.com/AmpliconSuite/CoRAL/wiki/Home/_edit#breakpoint-graphs).

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
concordant	chr7:55609190+->chr7:55609191-	2.697741	38
concordant	chr7:55610094+->chr7:55610095-	2.697741	41
concordant	chr7:56049369+->chr7:56049370-	4.150534	45
discordant	chr7:55610095-->chr7:55609190+	86.642611	869
discordant	chr7:56049369+->chr7:54763282-	85.189818	981
discordant	chr7:55155021-->chr7:55127266+	86.496697	978
...
PathConstraint: Path, Support
path_constraint e2+:1,c2-:1,e3+:1,c3-:1,e4+:1   6
path_constraint e4+:1,c4-:1,e5+:1,c5-:1,e6+:1   34
AmpliconIntervals: chr, start, end
interval        chr7    54659673        56149664
```
* Cycles file: 
For each amplicon, a tab-separated text file named ```output_dir_amplicon*_cycles.txt``` describing the list of cycles and paths returned from cycle extraction. Note that the cycles files output by CoRAL have mostly the same format as those output by [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (and therefore the files can be used interchangeably with AmpliconArchitect in most cases). Specifically a cycles file includes (i) the list of amplified intervals; (ii) the list of sequence edges; (iii) the list of cycles and paths, where an entry starts with ```0+``` and ends with ```0-``` in ```Segments``` indicates a path - these lines have the same format as AmpliconArchitect output. CoRAL's cycles files additionally include (iv) a list of longest (i.e., there are no paths that can form a sub/super-path to each other) path constraint indicated by long reads, and used in CoRAL's cycle extraction. Here is an example cycles file corresponding to the above graph file from GBM39.
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
Note that if ```--output-all-path-constraints``` is specified, then all path constraints given by long reads will be written to in ```*.cycles``` file.
* Other outputs include the ```output_dir_amplicon*_model.lp``` file(s) and ```output_dir_amplicon*_model.log``` file(s) given by Gurobi (integer program solver), for each amplicon, respectively describing the quadratic (constrained) program in a human readable format, and the standard output produced by Gurobi.


## 3. ```plot```
Usage: 
```coral plot <Required arguments> <Optional arguments>```

**3.1 Required arguments:**
If `--plot-graph` is given, `--graph` is required. If `--plot-cycles` is given `--cycles` is required.

| Argument                | Description                                                            |
|-------------------------|-----------------------------------------------------------------------|
| `--ref <choice>`        | Reference genome choice. Must be one of  `[hg19, hg38, GRCh38, mm10]` |
| `--bam <file>` | Bam file the run was based on                                         |
| `--graph <file>`        | AA-formatted `_graph.txt` file                                        |
| `--cycles <file>`       | AA-formatted `_cycles.txt` file                                       |
| `--output-dir <str>` | Directory for output files                                          |


**3.2 Optional arguments:**

| Argument                                   | Default                          | Description                                                                                                                               |
|--------------------------------------------|----------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| `--plot-graph`                             |                                  | Plot the AA graph file CN, SVs and coverage as a sashimi plot                                                                             |
| `--plot-cycles`                            |                                  | Plot the AA cycles file genome decompositions                                                                                             |
| `--only-cyclic-paths`                      |                                  | Only visualize the cyclic paths in the cycles file                                                                                        |
| `--num-cycles <int>`                       | `[all]`                          | Only plot the first `[arg]` cycles from the cycles file                                                                                   |
| `--max-coverage <float>`                   | `[1.25x max coverage in region]` | Do not extend coverage plot in graph sashimi plot above `[arg]` value                                                                     |
| `--min-mapq <int>`                         | 15                               | Do not use alignment in coverage plot with MAPQ value below `[arg]`                                                                       |
| `--gene-subset-list <str> <str> <str> ...` | `[all]`                          | Only indicate positions of the gene names in this list                                                                                    |
| `--hide-genes`                             |                                  | Do not plot positions of genes                                                                                                            |
| `--gene-fontsize <float>`                  | 12                               | Adjust fontsize of gene names                                                                                                             
| `--bushman-genes`                          |                                  | Only plot genes found in the [Bushman lab cancer-related gene list](http://www.bushmanlab.org/links/genelists) ('Bushman group allOnco'). | 
| `--region <chrom:pos1-pos2>`                | `[entire amplicon]`                | Only plot genome region in the interval given by `chrom:start-end`                                                                         |


## 4. ```hsr```
Usage: 
```coral hsr <Required arguments> <Optional arguments>```

**4.1 Required arguments:**

| Argument           | Descripion                                        |
|--------------------|---------------------------------------------------|
| `--lr-bam <file>`  | Coordinate-sorted and indexed long read .bam file |
| `--cycles <file>`  | AA-formatted `_cycles.txt` file                   |
| `--cn-segs <file>` | Long read segmented whole genome CN calls (.bed or CNVkit .cns file).            |
| `--normal-cov <float>` | Estimated coverage of diploid genome regions      |

**4.2 Optional arguments:**

| Argument                     | Default | Description                                                        |
|------------------------------|---------|--------------------------------------------------------------------|
| --bp_match_cutoff <int>      | 100     | Breakpoint matching cutoff distance (bp)                           |
| --bp_match_cutoff_clustering | 2000    | Crude breakpoint matching cutoff distance (bp) for clustering | 


## 5. ```cycle2bed```
CoRAL provides an option to convert its cycles output in AmpliconArchitect format ```*_cycles.txt``` into ```*.bed``` format (similar to [Decoil](https://github.com/madagiurgiu25/decoil-pre)), which makes it easier for downstream analysis of these cycles.

Usage: 
```coral cycle2bed <Required arguments> <Optional arguments>```

**5.1 Required arguments:**
* ```--cycle-file <FILE>``` - Input cycles file in AmpliconArchitect format.
* ```--output-file <FILE>```  - Output cycles file in ```*.bed``` format.

**5.2 Optional arguments:** 
* ```--num-cycles <INT>``` - If specified, only convert the first NUM_CYCLES cycles.

Here is an example output of ```cycle2bed``` given by the above cycles file from GBM39.
```
#chr	start	end	orientation	cycle_id	iscyclic	weight
chr7	54763282	55127266	+	1	True	82.346163
chr7	55155021	55609190	+	1	True	82.346163
chr7	55610095	56049369	+	1	True	82.346163
chr7	54763282	56049369	+	2	False	2.843655
```


## 6. ```cycle```
Usage: 
```coral cycle <Required arguments> <Optional arguments>```

**4.1 Required arguments:**

| Argument           | Descripion                                        |
|--------------------|---------------------------------------------------|
| `--graph <file>`   | AA-formatted `_graph.txt` file                   |
| `--output-dir <file>`  | Directory for output files                   |

**4.2 Optional arguments:**

| Argument                     | Default | Description                                                        |
|------------------------------|---------|--------------------------------------------------------------------|
| `--alpha <float>`      | 0.01     |  Parameter used to balance CN weight and path constraints in the objective function of greedy cycle extraction. Default value is 0.01, higher values favor the satisfaction of more path constraints.                           |
| `--solver-time-limit <int>` | 7200    | Time limit for cycle extraction (in seconds) | 
| `--threads <int>` | -1    | Number of threads for cycle extraction. If not specified, use all available cores. |
| `--solver <choice>` | gurobi_direct   | Solver for cycle extraction. Must be one of `[gurobi_direct, scip]` |
| `--output-all-path-constraints` | False    | If specified, output all path constraints given by long reads in `*_cycles.txt` file (see "Expected output" below). |
| `--postprocess-greedy-sol` | False    | If specified, automatically postprocess the cycles/paths returned in greedy cycle extraction, by solving the full quadratic program to minimize the number of cycles/paths starting with the greedy cycle extraction solution (as an initial solution). |


## FAQs
- `call_cnvs.sh` didn't produce segmented CN calls in a .cns file?
   - `cnvkit.py batch` contains multiple steps detailed in their 
   [documentation](https://cnvkit.readthedocs.io/en/stable/pipeline.html). The 
   errors from a particular stage don't always percolate up when running the
   complete pipeline via `batch`, so try running each stage separately to 
   pinpoint the root cause.
