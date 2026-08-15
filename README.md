#  CoRAL - Complete Reconstruction of Amplifications with Long reads
## Reference
CoRAL is a tool which utilizes aligned, single-molecule long-read data (.bam) as input, and identifies candidate ecDNA structures. The original Genome Research '24 paper is available here: https://genome.cshlp.org/content/34/9/1344.

**CoRAL only works on long-read whole-genome sequencing data (PacBio, Oxford Nanopore, etc.) - not targeted sequencing!**

## Installation
CoRAL can be installed and run on most modern Unix-like operating systems (e.g. Ubuntu 18.04+, CentOS 7+, macOS). Python >= 3.10 is required. Python 3.12 is recommended to reduce compatibility issues with some scientific dependencies.

1. Clone the repository
    ```bash
    git clone https://github.com/AmpliconSuite/CoRAL
    cd CoRAL
    ```

2. Install Poetry (if not already installed)
   ```bash
   pip install --user pipx
   pipx install poetry
   ```
   `pipx` installs Poetry in its own isolated environment, preventing conflicts with system Python packages. Do not use `pip install poetry` directly on a system Python.

   On macOS with a Homebrew-managed Python, `pip install --user pipx` will be rejected as an externally-managed environment (PEP 668). Use Homebrew instead:
   ```bash
   brew install pipx
   pipx ensurepath
   ```
   Either way, if `poetry` (or `pipx`-installed tools generally) isn't found on `PATH` afterward, run `pipx ensurepath` and restart your terminal (or `source` your shell config, e.g. `~/.zshrc` or `~/.bashrc`) to pick up the updated `PATH`.

3. Install CoRAL dependencies
   ```bash
   poetry install
   ```
   Poetry creates an isolated virtual environment automatically. If `pysam` fails to build, install `htslib` first (`sudo apt install libhtslib-dev` on Debian/Ubuntu, or `brew install htslib` on macOS), then re-run `poetry install`.

   `cvxopt` may fall back to a source build that fails with `fatal error: 'umfpack.h' file not found`. This means it needs SuiteSparse, which provides UMFPACK. Install it first (`sudo apt install libsuitesparse-dev` on Debian/Ubuntu, or `brew install suite-sparse` on macOS), then re-run `poetry install`.

4. Activate the environment and verify the installation
   ```bash
   source $(poetry env info --path)/bin/activate
   coral --help
   ```
   This works with all versions of Poetry. To deactivate the environment, run `deactivate`.

### Re-activating the environment
After the initial install, you only need to re-run the activation command each time you start a new shell session (e.g. after logging out and back in):
```bash
cd /your/path/to/CoRAL/
source $(poetry env info --path)/bin/activate
```
Run this from the CoRAL repository directory. To avoid doing this manually each session, you can add it to your `~/.bashrc` or `~/.zshrc`.

5. [Download a Gurobi optimizer license](https://support.gurobi.com/hc/en-us/articles/360040541251-How-do-I-obtain-a-free-academic-license) (free for academic use)
   - Place the `gurobi.lic` file in `$HOME/gurobi.lic`.

6. Finish installing CNVkit dependencies (recommended)
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


- If you don't have these then you can run CNVkit (installed as a dependency) to generate them, by running

   `./scripts/call_cnvs.sh <input.bam> ./reference/hg38full_ref_5k.cnn <output_dir>`

   This will create a file called `[input].cns`, which you can feed to CoRAL for it's `--cn_segs` argument.

## Running CoRAL

CoRAL is run as a series of **stages**, each invoked as a separate mode:

`coral [mode] [mode arguments]`

A standard analysis is **stage 1 → stage 2**, since `reconstruct` requires the seed intervals produced by `seed`. Stages 3-5 are optional and operate on the outputs of stage 2.

| Stage | Mode | Purpose |
|-------|------|---------|
| 1 | [`seed`](#1-seed) | Identify seed amplified intervals (copy number gain regions) from genome-wide CN calls. **Run first** - its output is a required input to `reconstruct`. |
| 2 | [`reconstruct`](#2-reconstruct) | **The main running mode of CoRAL.** Build the breakpoint graph from the BAM and extract cycles, producing `*_graph.txt` and `*_cycles.txt`. |
| 3 | [`cycle`](#3-cycle) | *Optional.* Re-run **only** cycle extraction on an existing breakpoint graph, without re-reading the BAM. |
| 4 | [`plot`](#4-plot) | *Optional.* Visualize the reconstructed cycles/paths and/or the breakpoint graph. |
| 5 | [`hsr`](#5-hsr) | *Optional, experimental.* Identify candidate chromosomal homogeneously staining region (HSR) integration points for ecDNA. |

Once reconstruction is complete, the resulting amplicons can be classified by amplification type (ecDNA, BFB, etc.) with AmpliconClassifier - see [Amplicon classification](#6-amplicon-classification). Format-conversion helpers are described under [Optional utilities](#optional-utilities).


## 1. ```seed```
As the seed amplification intervals are required by the main script ```reconstruct``` mode, it is suggested the user first run ```seed``` mode to generate seed amplification intervals.

Usage: 
```coral seed <Required arguments> <Optional arguments>```

**Required arguments:**
* ```--cn-seg <file>``` - Long read segmented whole genome CN calls. Accepts CNVkit `.cns` format, or a tab-separated `.bed` file where the **last column** is the copy number value.
* ```--output-prefix <string>``` - Prefix of the output ```*_CNV_SEEDS.bed``` file.
* ```--lr-bam <file>``` - Coordinate sorted BAM file, used to read chromosome lengths from the header.
* ```--centromere-file <file>``` - Centromere BED file for the reference genome used.

  The file must be tab-separated with columns `chr`, `start`, `end` (BED3 format; additional columns are ignored). Each contig may have one entry, or two entries that overlap or directly abut (they will be merged). Multiple distinct non-adjacent regions for the same contig will raise an error. Contigs absent from the file are treated as having no centromere (all segments on that contig are eligible as seeds). Example:
  ```
  chr1	121700000	125100000
  chr2	91800000	100200000
  ```
  The UCSC `centromeres` or `gap` track (available via the Table Browser for most assemblies) is a convenient source for this file.

**Optional arguments:**
* ```--gain <float>``` - A minimum CN threshold (with the assumption of diploid genome) for a particular CN segment to be considered as a seed. Default is 6.0.
* ```--min-seed-size <int>``` - Minimum size (in bp) for a CN segment to be considered as a seed. Default is 100000.
* ```--max-seg-gap <int>``` - Maximum gap size (in bp) to merge two proximal CN segments to be considered as seed intervals. If at least two segments are merged, then they will be treated as a single candidate to be filtered with ```--min-seed-size```, and their aggregate size will be compared with the value. Default is 300000.
* ```--extra-contigs <file>``` - Plain-text file of additional contig names (one per line) to include alongside standard chromosomes when reading chromosome sizes from the BAM header.


## 2. ```reconstruct```
Usage: 
```reconstruct <Required arguments> <Optional arguments>```

**2.1 Required arguments:**
* ```--lr-bam <file>``` - Coordinate sorted ```*.BAM``` file, with ```*.bai``` index (mapped to the provided reference genome) in the same directory.
* ```--cnv-seed <file>``` - ```*.bed``` file with a putative list of seed amplification intervals. The seed amplification intervals can be obtained through [running ```seed``` mode](#CoRAL.py-```seed```), or provided manually.
* ```--output-prefix <path>``` - Prefix (including directory) to which the output ```graph.txt``` and ```cycles.txt``` files will be written.
* ```--cn-seg <file>``` - Long read segmented whole genome CN calls. Accepts CNVkit `.cns` format, or a tab-separated `.bed` file where the **last column** is the copy number value.

**2.2 Optional arguments:**
* ```--min-bp-support <float>``` - Filter out breakpoints with less than (min_bp_support * normal coverage) long read support in breakpoint graph construction. The default value is set to 1.0, meaning to filter out breakpoints supported by less than diploid coverage, but ***it is highly recommended to specify a much larger value, e.g. 10.0 to obtain a cleaner breakpoint graph and the dominating ecDNA cycle(s).***
* ```--skip-cycle-decomp``` - If specified, will stop by only outputting the breakpoint graph files ```*_graph.txt``` (see [**Expected output**](#2.3-Expected-output) below) for all amplicons and not extract cycles from the graph and output ```*_cycles.txt```.
* ```--output-path_constraints <all|longest|none>``` - Options to output path constraints given by long reads. By default, CoRAL will only output the longest path constraints in ```*_graph.txt``` and ```*_cycles.txt``` files (see "Expected output" below). If ```all``` or ```none``` is specified, output all or no path constraints in both graph and cycles files.
* ```--cycle-decomp-mode <min_cycles|max_weight> ``` - min_cycles: minimize the number of cycles/paths, if solver could not find a feasible solution within time limit, switch to greedy cycle extraction; max_weight: greedy cycle extraction by iteratively extracting cycles or paths with maximum length weighted CN and satisfying the most path constraints. Default mode is ```max_weight```.
* ```--cycle-decomp-alpha <float between [0, 1]>``` - Parameter used to balance CN weight and path constraints in the objective function of greedy cycle extraction. Default value is 0.01, higher values favor the satisfaction of more path constraints.
* ```--solver-time-limit <int>``` - Maximum running time (in seconds) reserved for solving a single quadratic program using the chosen integer program solver (e.g. Gurobi, SCIP). The solver would return the best solution(s) it currently found, regardless of the optimality status, when reaching this time limit. Default value is 7200 (i.e., 2 hours).
* ```--solver-threads <int>``` - Number of threads reserved for for solving the quadratic program with Gurobi (integer program solver). If not specified (and by default), the solver would attempt to use up all available cores in the working machine. 
* ```--solver <choice>``` - Solver for cycle extraction. Must be one of `[gurobi_direct, scip]`.
* ```--global-time-limit <int>``` - Maximum running time (in seconds) reserved for the entire cycle extraction process. Default value is 21600 (i.e., 6 hours).
* ```--postprocess-greedy-sol``` - If specified, automatically postprocess the cycles/paths returned in greedy cycle extraction, by solving the full quadratic program to minimize the number of cycles/paths starting with the greedy cycle extraction solution (as an initial solution).
*	```--log-file <file>``` - Name of the main ```*.log``` file, which can be used to trace the status of ```reconstruct``` run(s).
* ```--verbose```/```-v``` - Write full debug output to the ```*.log``` file, instead of only ```INFO```-level messages and above. This typically increases the size of the log by 10-50x, and is mainly useful for gathering detail when reporting an issue.
* ```--cache-reads``` - Write a ```<output-prefix>_chimeric_alignments.pickle``` cache of the chimeric reads fetched from the ```*.BAM```, **which can exceed 1 GB on high-coverage samples**. If that cache file is already present, a subsequent run using the same ```--output-prefix``` will load it and skip re-scanning the BAM, which can substantially reduce runtime. The cache records which BAM it was built from and is ignored automatically if that BAM changes.

**2.3 Expected output:**

CoRAL may identify and reconstruct several distinct focal amplifications in the input ```*.BAM``` sample. Each is organized as an *amplicon* - a connected component of amplified intervals and their connections by discordant edges - and is numbered `amplicon1`, `amplicon2`, and so on.

All output files are named using the value given to ```--output-prefix```, which includes the output directory (e.g. `--output-prefix results/GBM39` writes `results/GBM39_amplicon1_graph.txt`). The files written are:

| File | Written per | Contents |
|------|-------------|----------|
| `<prefix>_amplicon<N>_graph.txt` | amplicon | **The breakpoint graph.** Sequence, concordant, and discordant edges with predicted copy number, plus path constraints and amplicon intervals. |
| `<prefix>_amplicon<N>_cycles.txt` | amplicon | **The reconstructed cycles and paths** with their copy counts, and which path constraints each satisfies. Not written if ```--skip-cycle-decomp``` is given. |
| `<prefix>_summary.txt` | run | Per-amplicon summary of the reconstruction (intervals, amplicon count, whether cycle extraction succeeded). Used by downstream tools to locate the run. |
| `<prefix>_reconstruct.log` | run | Run log. See ```--log-file``` and ```--verbose``` above. |
| `<prefix>_amplicon<N>_model.lp` / `.log` | amplicon | The integer program handed to the solver, in human-readable form, and the solver's own output. Diagnostic only. |

The two files that matter for downstream analysis are `*_graph.txt` and `*_cycles.txt`. Both use the [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (AA) format and are interchangeable with AA's equivalents, which is what allows them to be passed directly to [`plot`](#4-plot), [`cycle`](#3-cycle), and [AmpliconClassifier](#6-amplicon-classification).

The remainder of this section describes the two formats in detail.

* Graph file: For each amplicon, a tab-separated text file named ```<prefix>_amplicon*_graph.txt``` describing the *sequence edges*, *concordant edges* and *discordant edges* in the graph and their predicted copy count. Note that the graph files outputted by CoRAL have the same format as those outputted by [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (and therefore the files can be used interchangeably with AmpliconArchitect). Here is an example graph file from GBM39, a cell line with *EGFR* amplified on ecDNA.
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
For each amplicon, a tab-separated text file named ```<prefix>_amplicon*_cycles.txt``` describing the list of cycles and paths returned from cycle extraction. Note that the cycles files output by CoRAL have mostly the same format as those output by [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (and therefore the files can be used interchangeably with AmpliconArchitect in most cases). Specifically a cycles file includes (i) the list of amplified intervals; (ii) the list of sequence edges; (iii) the list of cycles and paths, where an entry starts with ```0+``` and ends with ```0-``` in ```Segments``` indicates a path - these lines have the same format as AmpliconArchitect output. CoRAL's cycles files additionally include (iv) a list of longest (i.e., there are no paths that can form a sub/super-path to each other) path constraint indicated by long reads, and used in CoRAL's cycle extraction. Here is an example cycles file corresponding to the above graph file from GBM39.
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


## 3. ```cycle```

`reconstruct` already performs cycle extraction, so **most users never need this mode.** It exists to re-run *only* the cycle extraction half of `reconstruct`, starting from a breakpoint graph that already exists. Because it reads the `*_graph.txt` file rather than the BAM, it skips the expensive alignment-scanning step entirely. Use it when you want to:

- re-extract cycles under different parameters (e.g. a different ```--cycle-decomp-mode```, ```--alpha```, or a longer ```--solver-time-limit```) without rebuilding the graph;
- resume a run that used ```--skip-cycle-decomp```;
- extract cycles from a graph produced by [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) rather than by CoRAL.

**NOTE:** a CoRAL-generated graph must come from CoRAL v2.1.0 or later, since the `path constraints` and `amplicon intervals` sections that this mode requires were only added to `*_graph.txt` in that release. Graphs from earlier versions must be regenerated with `reconstruct`.

Usage: 
- ```coral cycle <Required arguments> <Optional arguments>``` for cycle extraction from a single amplicon (breakpoint graph);
- ```coral cycle_all <Required arguments> <Optional arguments>``` for cycle extraction from all amplicons in a directory.

**3.1 Required arguments:**

| Argument           | Descripion                                        |
|--------------------|---------------------------------------------------|
| `--graph <file>` (for `cycle` mode)  | AA or CoRAL-formatted `_graph.txt` file.          |
| `--bp-dir <directory>`  (for `cycle_all` mode)  | Directory containing AA or CoRAL-formatted `_graph.txt` files for different amplicons. |
| `--output-prefix <path>`  | Prefix (including directory) of the output ```cycles.txt``` files.                    |

**3.2 Optional arguments:**

| Argument                     | Default | Description                                                        |
|------------------------------|---------|--------------------------------------------------------------------|
| `--alpha <float>`      | 0.01     |  Parameter used to balance CN weight and path constraints in the objective function of greedy cycle extraction. Default value is 0.01, higher values favor the satisfaction of more path constraints.                           |
| `--solver-time-limit <int>` | 7200    | Time limit for cycle extraction (in seconds) | 
| `--threads <int>` | -1    | Number of threads for cycle extraction. If not specified, use all available cores. |
| `--solver <choice>` | gurobi_direct   | Solver for cycle extraction. Must be one of `[gurobi_direct, scip]`. |
| `--cycle-decomp-mode <choice>` | max_weight | Cycle extraction mode. Must be one of `[min_cycles, max_weight]` (see `reconstruct` mode). |
| `--output-path-constraints <choice>` | longest | Options for output path constraints in graph and cycle files. Must be one of `[all, longest, none]` (see `reconstruct` mode). |
| `--postprocess-greedy-sol` | False    | If specified, automatically postprocess the cycles/paths returned in greedy cycle extraction, by solving the full quadratic program to minimize the number of cycles/paths starting with the greedy cycle extraction solution (as an initial solution). |


## 4. ```plot```
Usage: 
```coral plot <Required arguments> <Optional arguments>```

**4.1 Required arguments:**

| Argument                | Description                                                            |
|-------------------------|-----------------------------------------------------------------------|
| `--ref <choice>`        | Reference genome. Must be one of `[hg19, hg38, mm10]`                 |
| `--output-prefix <str>` | Prefix of output files                                                 |

At least one of `--graph` or `--cycles` must also be provided.

**4.2 Optional arguments:**

| Argument                                   | Default                          | Description                                                                                                                               |
|--------------------------------------------|----------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| `--graph <file>`                           |                                  | AA-formatted `_graph.txt` file                                                                                                            |
| `--cycles <file>`                          |                                  | AA-formatted `_cycles.txt` file                                                                                                           |
| `--bam <file>`                             |                                  | Sorted and indexed BAM file. If provided, the coverage track (gray bars) is computed from the aligned reads in small windows and show depth variation within each segment; the track uses coverage stored in the graph file, one value per segment. |
| `--only-cyclic-paths`                      |                                  | Only visualize the cyclic paths in the cycles file                                                                                        |
| `--num-cycles <int>`                       | `[all]`                          | Only plot the first `[arg]` cycles from the cycles file                                                                                   |
| `--max-coverage <float>`                   | `[1.25x max coverage in region]` | Do not extend coverage plot in graph sashimi plot above `[arg]` value                                                                     |
| `--min-mapq <int>`                         | 15                               | Do not use alignment in coverage plot with MAPQ value below `[arg]`                                                                       |
| `--gene-subset-list <str> <str> <str> ...` | `[all]`                          | Only indicate positions of the gene names in this list                                                                                    |
| `--gene-subset-file <file>`                | `[all]`                          | Text or CSV file of gene names to indicate in the plot. Gene names may be newline, whitespace, or comma delimited. Empty files emit a warning and fall back to plotting all genes unless `--gene-subset-list` is also provided. |
| `--hide-genes`                             |                                  | Do not plot positions of genes                                                                                                            |
| `--font-size <float>`                      | 1                                | Multiply all plot text and axis styling by this non-negative value. This includes titles, axis labels, tick labels, tick marks, axis lines, legends, cycle labels, and gene names. A value of `0` hides all text and axis ticks while retaining plotted data and gene tracks. When explicitly supplied, this option overrides `--gene-fontsize`. |
| `--gene-fontsize <float>`                  | 12                               | Set the gene-name font size when `--font-size` is not explicitly supplied. |
| `--bushman-genes`                          |                                  | Only plot genes found in the [Bushman lab cancer-related gene list](http://www.bushmanlab.org/links/genelists) ('Bushman group allOnco'). | 
| `--region <chrom:pos1-pos2>`                | `[entire amplicon]`                | Only plot genome region in the interval given by `chrom:start-end`                                                                         |

Graph plots can be generated from a graph file alone:

```bash
coral plot --ref hg38 --graph results/GBM39_amplicon1_graph.txt --output-prefix sample_graph
```

In this mode, CoRAL uses the predicted CN and average coverage columns already present in `_graph.txt`. If a BAM is provided, CoRAL preserves the existing behavior of extracting the gray coverage bars directly from alignments:

```bash
coral plot --ref hg38 --graph results/GBM39_amplicon1_graph.txt --bam sample.bam --output-prefix sample_graph_bam
```

Graph plots write a separate legend figure next to the requested output prefix. The legend documents coverage source, predicted segment CN, discordant-edge orientation colors, and that discordant-edge width scales with discordant read count.

Gene subsets may be supplied directly or through a file:

```bash
coral plot --ref hg38 --graph results/GBM39_amplicon1_graph.txt --gene-subset-file genes.csv --output-prefix sample_graph_genes
```

All plot text and axis styling can be scaled together:

```bash
# Twice the default text, label, tick, and axis sizes.
coral plot --ref hg38 --graph sample_data/test4/amplicon1_graph.txt --font-size 2 --output-prefix sample_graph_large_text

# Half the default text, label, tick, and axis sizes.
coral plot --ref hg38 --graph sample_data/test4/amplicon1_graph.txt --font-size 0.5 --output-prefix sample_graph_small_text

# Keep plotted data and gene tracks, but omit all text and axis ticks.
coral plot --ref hg38 --graph sample_data/test4/amplicon1_graph.txt --font-size 0 --output-prefix sample_graph_no_text
```

`--font-size` is also accepted by `coral plot_all` and applies to every generated plot.

If both `--font-size` and `--gene-fontsize` are explicitly supplied,
`--font-size` takes precedence and CoRAL warns that `--gene-fontsize` is being
ignored. Omitting `--font-size` preserves the existing `--gene-fontsize`
behavior.


## 5. ```hsr```

> **Experimental.** Treat the reported integration points as candidates for follow-up rather than as confident calls.

Identifies candidate locations where an ecDNA may have integrated into a chromosome as a homogeneously staining region (HSR).

Usage: 
```coral hsr <Required arguments> <Optional arguments>```

**5.1 Required arguments:**

| Argument           | Descripion                                        |
|--------------------|---------------------------------------------------|
| `--lr-bam <file>`  | Coordinate-sorted and indexed long read .bam file |
| `--cycles <file>`  | AA-formatted `_cycles.txt` file                   |
| `--cn-seg <file>` | Long read segmented whole genome CN calls. Accepts CNVkit `.cns` format, or a tab-separated `.bed` file where the **last column** is the copy number value.            |
| `--normal-cov <float>` | Estimated coverage of diploid genome regions      |

**5.2 Optional arguments:**

| Argument                          | Default | Description                                                        |
|-----------------------------------|---------|--------------------------------------------------------------------|
| `--bp-match-cutoff <int>`         | 100     | Breakpoint matching cutoff distance (bp)                           |
| `--bp-match-cutoff-clustering <int>` | 2000 | Crude breakpoint matching cutoff distance (bp) for clustering      |
| `--extra-contigs <file>`          |         | Plain-text file of additional contig names (one per line) to include alongside standard chromosomes when reading chromosome sizes from the BAM header. |


## 6. Amplicon classification

CoRAL reconstructs amplicon structures but does not label them by mechanism. To determine whether a reconstructed amplicon is **ecDNA**, **BFB**, **complex non-cyclic**, or **linear**, pass its `*_graph.txt` and `*_cycles.txt` files to [AmpliconClassifier](https://github.com/AmpliconSuite/AmpliconClassifier) (AC).

**CoRAL outputs are supported by AC v2.0.0 and later.** Earlier AC releases predate CoRAL support and should not be used.

**6.1 Installing AmpliconClassifier**

AC is a separate tool with its own environment:
```bash
conda create -n ampliconclassifier "python>=3.8"
conda activate ampliconclassifier
git clone https://github.com/AmpliconSuite/AmpliconClassifier.git
cd AmpliconClassifier
python -m pip install -e .
```

AC also needs the AA data repository, which supplies its reference annotations (gene models, oncogene list, mappability and centromere tracks):
```bash
echo 'export AA_DATA_REPO=~/data_repo/' >> ~/.bashrc
source ~/.bashrc
mkdir -p $AA_DATA_REPO

ref=GRCh38   # or hg19, GRCh37, mm10, GRCh38_viral
wget https://refs.ampliconrepository.org/data/module_support_files/AmpliconArchitect/${ref}.tar.gz -P $AA_DATA_REPO
tar -xzf $AA_DATA_REPO/${ref}.tar.gz -C $AA_DATA_REPO
rm $AA_DATA_REPO/${ref}.tar.gz
```
Use the plain build name (`GRCh38`), not the `_indexed` variant - the latter adds a BWA index that AC never uses. See the [AC README](https://github.com/AmpliconSuite/AmpliconClassifier#1-installation) for the full instructions, including optional dependencies.

**6.2 Classifying CoRAL output**

Point AC at the directory containing your CoRAL results; it will find the graph/cycles pairs itself and classify every amplicon:
```bash
conda activate ampliconclassifier
amplicon_classifier.py --ref GRCh38 --AA_results /path/to/coral/output/ -o my_run > my_run_classifier.log
```

A few notes specific to CoRAL:
* `--ref` takes an **AA** build name. A CoRAL run against `hg38` corresponds to `--ref GRCh38` (both use `chr`-prefixed names); `hg19` and `mm10` are passed through unchanged.
* AC pairs each `*_cycles.txt` with the `*_graph.txt` sharing the same prefix, so keep both files together and leave them under the names CoRAL wrote. It also picks up CoRAL's `*_summary.txt`.
* To classify a single amplicon instead of a whole directory, use `--cycles`/`--graph`:
  ```bash
  amplicon_classifier.py --ref GRCh38 \
      --cycles results/GBM39_amplicon1_cycles.txt \
      --graph  results/GBM39_amplicon1_graph.txt \
      -o GBM39_amplicon1
  ```

The primary result is `my_run_amplicon_classification_profiles.tsv`, giving each amplicon an `amplicon_decomposition_class` along with `ecDNA+`, `BFB+`, and `FAN+` flags. AC also writes a gene list, per-feature property and complexity tables, and a summary results table. See the [AC README](https://github.com/AmpliconSuite/AmpliconClassifier#3-outputs) for the full description of each output.

**6.3 Sharing results on AmpliconRepository**

Classified results can be deposited in [AmpliconRepository.org](https://ampliconrepository.org), a public portal for browsing and sharing focal amplification calls. AC output is **required** - the repository uses AC's `*_result_table.tsv` to decide which samples to ingest, and recognizes CoRAL runs by their `*_summary.txt` files.

Package the CoRAL and AC outputs together as a `.tar.gz` or `.zip`, excluding sequence and working files:
```bash
tar -czf my_study.tar.gz \
    --exclude='*.bam*' --exclude='*.cram*' --exclude='*.fastq*' --exclude='*.fq*' \
    --exclude='*.pickle' --exclude='*.cnn' --exclude='*.cnr' \
    my_study/
```
Then log in, choose **New Project**, and upload. Full instructions are in the [AmpliconRepository getting-started guide](https://docs.ampliconrepository.org/en/latest/getting-started/).


## Optional utilities

These helpers are not part of the reconstruction pipeline and can be run at any time on existing output.

### ```cycle2bed```
Converts a cycles file in AmpliconArchitect format ```*_cycles.txt``` into ```*.bed``` format (similar to [Decoil](https://github.com/madagiurgiu25/decoil-pre)), which makes it easier for downstream analysis of these cycles.

Usage: 
```coral cycle2bed <Required arguments> <Optional arguments>```

**Required arguments:**
* ```--cycle-file <file>``` - Input cycles file in AmpliconArchitect format.
* ```--output-file <file>```  - Output cycles file in ```*.bed``` format.

**Optional arguments:** 
* ```--num-cycles <int>``` - If specified, only convert the first NUM_CYCLES cycles.

Here is an example output of ```cycle2bed``` given by the above cycles file from GBM39.
```
#chr	start	end	orientation	cycle_id	iscyclic	weight
chr7	54763282	55127266	+	1	True	82.346163
chr7	55155021	55609190	+	1	True	82.346163
chr7	55610095	56049369	+	1	True	82.346163
chr7	54763282	56049369	+	2	False	2.843655
```

## FAQs
- `call_cnvs.sh` didn't produce segmented CN calls in a .cns file?
   - `cnvkit.py batch` contains multiple steps detailed in their 
   [documentation](https://cnvkit.readthedocs.io/en/stable/pipeline.html). The 
   errors from a particular stage don't always percolate up when running the
   complete pipeline via `batch`, so try running each stage separately to 
   pinpoint the root cause.
