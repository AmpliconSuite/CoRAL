import re

COMMA_NUMBER_REGEX = r"\d{1,3}(?:,\d{3})*"
AMPLICON_SEPARATOR = "-----------------------------------------------"


## Summary shared patterns

# Header patterns
VERSION_TEMPLATE = "CoRAL v{version}"
VERSION_PATTERN = re.compile(r"CoRAL v(\d+\.\d+\.\d+)")
PROFILE_ENABLED_TEMPLATE = "Profiling Enabled: {enabled}"
PROFILE_ENABLED_PATTERN = re.compile(r"Profiling Enabled: (True|False)")


# Solver patterns
SOLVER_PATTERN = re.compile(r"Solver: (\S+)")
THREADS_PATTERN = re.compile(r"Threads: (\d+)")
TIME_LIMIT_PATTERN = re.compile(r"Time Limit: (\d+) s")


# Per-amplicon summary patterns
AMPLICON_SIZE_PATTERN = re.compile(f"Amplicon Size: ({COMMA_NUMBER_REGEX})")
AMPLICON_BREAKPOINTS_PATTERN = re.compile(
    f"# Discordant Edges: ({COMMA_NUMBER_REGEX})"
)

CYCLE_DECOMP_STATUS_TEMPLATE = "Cycle Decomposition Status: {status}"
CYCLE_DECOMP_STATUS_PATTERN = re.compile(r"Cycle Decomposition Status: (\S+)")
