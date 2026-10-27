"""Tests for resolving the running CoRAL version.

Output headers must report the version of the code that actually ran. When a
source checkout is imported (PYTHONPATH or editable install) while a different
CoRAL release is installed in the environment, the checkout's pyproject.toml
wins; otherwise the installed distribution metadata is used.
"""

from __future__ import annotations

import importlib.metadata
import pathlib
from collections.abc import Iterator

import pytest

from coral import version as coral_version


@pytest.fixture(autouse=True)
def clear_version_cache() -> Iterator[None]:
    coral_version.get_coral_version.cache_clear()
    yield
    coral_version.get_coral_version.cache_clear()


def _fake_package(tmp_path: pathlib.Path, pyproject_text: str | None) -> pathlib.Path:
    pkg = tmp_path / "coral"
    pkg.mkdir()
    module_file = pkg / "version.py"
    module_file.write_text("")
    if pyproject_text is not None:
        (tmp_path / "pyproject.toml").write_text(pyproject_text)
    return module_file


def test_source_tree_version_matches_repo_pyproject() -> None:
    repo_pyproject = pathlib.Path(coral_version.__file__).resolve().parent.parent / "pyproject.toml"
    if not repo_pyproject.is_file():
        pytest.skip("not running from a source checkout")
    assert coral_version.get_coral_version() in repo_pyproject.read_text()


def test_source_tree_preferred_over_installed_metadata(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module_file = _fake_package(
        tmp_path, '[tool.poetry]\nname = "CoRAL"\nversion = "9.9.9"\n'
    )
    monkeypatch.setattr(coral_version, "__file__", str(module_file))
    monkeypatch.setattr(importlib.metadata, "version", lambda _: "1.0.0")
    assert coral_version.get_coral_version() == "9.9.9"


def test_falls_back_to_installed_metadata_without_pyproject(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module_file = _fake_package(tmp_path, None)
    monkeypatch.setattr(coral_version, "__file__", str(module_file))
    monkeypatch.setattr(importlib.metadata, "version", lambda _: "1.0.0")
    assert coral_version.get_coral_version() == "1.0.0"


def test_ignores_pyproject_of_another_project(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module_file = _fake_package(
        tmp_path, '[tool.poetry]\nname = "something-else"\nversion = "9.9.9"\n'
    )
    monkeypatch.setattr(coral_version, "__file__", str(module_file))
    monkeypatch.setattr(importlib.metadata, "version", lambda _: "1.0.0")
    assert coral_version.get_coral_version() == "1.0.0"


def test_unknown_when_nothing_resolves(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module_file = _fake_package(tmp_path, None)
    monkeypatch.setattr(coral_version, "__file__", str(module_file))

    def _missing(_: str) -> str:
        raise importlib.metadata.PackageNotFoundError("coral")

    monkeypatch.setattr(importlib.metadata, "version", _missing)
    assert coral_version.get_coral_version() == "unknown"
