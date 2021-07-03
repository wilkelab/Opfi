# Taken from the pytest documentaion site:
# https://docs.pytest.org/en/latest/example/

import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )
    parser.addoption(
        "--runprop", action="store_true", default=False, help="run very slow property-based tests"
    )
    parser.addoption(
        "--runmmseqs", action="store_true", default=False, help="run tests that have mmseqs2 as a dependency"
    )
    parser.addoption(
        "--rundiamond", action="store_true", default=False, help="run tests that have diamond as a dependency"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")
    config.addinivalue_line("markers", "proptest: mark test as a very slow property test")
    config.addinivalue_line("markers", "mmseqs: mark test as having mmseqs2 as a dependency")
    config.addinivalue_line("markers", "diamond: mark test as having diamond as a dependency")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        skip_slow = pytest.mark.skip(reason="need --runslow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)
    if not config.getoption("--runprop"):
        skip_prop = pytest.mark.skip(reason="need --runprop option to run")
        for item in items:
            if "proptest" in item.keywords:
                item.add_marker(skip_prop)
    if not config.getoption("--runmmseqs"):
        skip_mmseqs = pytest.mark.skip(reason="need --runmmseqs option to run")
        for item in items:
            if "mmseqs" in item.keywords:
                item.add_marker(skip_mmseqs)
    if not config.getoption("--rundiamond"):
        skip_diamond = pytest.mark.skip(reason="need --rundiamond option to run")
        for item in items:
            if "diamond" in item.keywords:
                item.add_marker(skip_diamond)
