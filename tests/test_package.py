"""Package and optional-integration tests."""

import subprocess
import sys
from types import SimpleNamespace

import pytest

import md_davis
from md_davis.integrations.pymol import (
    PyMOLUnavailableError,
    initialize_pymol,
    require_pymol,
)


def test_package_import_is_lightweight():
    script = (
        "import sys; import md_davis; "
        "assert 'md_davis.gui' not in sys.modules; "
        "assert 'pymol' not in sys.modules; "
        "print(md_davis.__version__)"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip() == md_davis.__version__ == "0.4.1"


def test_pymol_error_explains_open_source_install(monkeypatch):
    monkeypatch.setitem(sys.modules, "pymol", None)
    with pytest.raises(PyMOLUnavailableError, match="pymol-open-source"):
        require_pymol()


@pytest.mark.parametrize(
    ("window", "expected_args"),
    [(False, ["pymol", "-cQ", "--test"]), (True, ["pymol", "--test"])],
)
def test_open_source_pymol_initialization(monkeypatch, window, expected_args):
    fake_pymol = SimpleNamespace(cmd=object(), launched=False)

    def finish_launching():
        fake_pymol.launched = True

    fake_pymol.finish_launching = finish_launching
    monkeypatch.setitem(sys.modules, "pymol", fake_pymol)

    command = initialize_pymol(window=window, extra_args=["--test"])

    assert command is fake_pymol.cmd
    assert fake_pymol.pymol_argv == expected_args
    assert fake_pymol.launched is True
