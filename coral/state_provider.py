"""Shared state used across various submodules - either injected from config or set on initialization."""

import time

# Time since main function starts. Used for logging output
TSTART: float = 0.0


def get_timed_task_log(task_desc: str) -> str:
    return f"#TIME {time.time() - TSTART:.4f}s: {task_desc}"
