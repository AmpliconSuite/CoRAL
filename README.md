#  CoRAL - Complete Reconstruction of Amplifications with Long reads
## Reference
https://www.biorxiv.org/content/10.1101/2024.02.15.580594v1

## Installation
CoRAL can be installed and run on most modern Unix-like operating systems (e.g. Ubuntu 18.04+, CentOS 7+, macOS). 

CoRAL requires python>=3.7, and the following packages.
* pysam(>=0.1.7) https://pysam.readthedocs.io/en/stable/ for reading mapped sequences in ```*.BAM``` format
* cvxopt https://cvxopt.org/ for estimating CN in breakpoint graph.
* Gurobi (>=9.1, for Python) https://www.gurobi.com/documentation/current/refman/py_python_api_overview.html for solving the quadratic (constrained) program to extract cycles/paths from the breakpoint graph.
* (Optional) CNVkit 

## Command line arguments to run CoRAL

