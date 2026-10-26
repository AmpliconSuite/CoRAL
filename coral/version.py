"""Resolve the version of the CoRAL code that is actually running."""

from __future__ import annotations

import functools
import importlib.metadata
import logging
import pathlib
import re

logger = logging.getLogger(__name__)

_PYPROJECT_VERSION_RE = re.compile(r'^\s*version\s*=\s*"([^"]+)"', re.MULTILINE)


def _version_from_source_tree() -> str | None:
    """Read the version from the pyproject.toml beside the imported package.

    Returns None when CoRAL is not being run from a source checkout (e.g. a
    regular site-packages install, where no pyproject.toml sits next to it).
    """
    pyproject = pathlib.Path(__file__).resolve().parent.parent / "pyproject.toml"
    if not pyproject.is_file():
        return None
    try:
        text = pyproject.read_text()
    except OSError:
        return None

    try:
        import tomllib

        poetry = tomllib.loads(text).get("tool", {}).get("poetry", {})
        if str(poetry.get("name", "")).lower() != "coral":
            return None
        version = poetry.get("version")
        return str(version) if version else None
    except ImportError:
        # Python < 3.11: no tomllib, fall back to a plain scan of the file.
        if not re.search(r'^\s*name\s*=\s*"coral"', text, re.MULTILINE | re.IGNORECASE):
            return None
        match = _PYPROJECT_VERSION_RE.search(text)
        return match.group(1) if match else None
    except Exception:
        return None


@functools.lru_cache(maxsize=1)
def get_coral_version() -> str:
    """Version of the running CoRAL code.

    Prefers the source tree the `coral` package was imported from, so a
    checkout run via PYTHONPATH (or an editable install) reports its own
    version rather than whatever distribution is installed in the environment.
    Falls back to the installed distribution metadata, then to "unknown".
    """
    version = _version_from_source_tree()
    if version:
        return version
    try:
        return importlib.metadata.version("coral")
    except importlib.metadata.PackageNotFoundError:
        logger.warning("Could not determine CoRAL version; reporting 'unknown'.")
        return "unknown"
