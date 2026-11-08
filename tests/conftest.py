"""Load the top-level scripts (one has a hyphen in its name) as modules."""
import importlib.util
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def sb():
    return _load("spot_builder", "spot-builder.py")


@pytest.fixture(scope="session")
def tar():
    return _load("taratura", "taratura.py")
