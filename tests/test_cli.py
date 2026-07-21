"""Tests for the argparse command-line interface."""

import argparse
import subprocess
import sys

from md_davis.cli import build_parser


def run_cli(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "md_davis", *arguments],
        capture_output=True,
        check=False,
        text=True,
    )


def test_cli_uses_argparse():
    assert isinstance(build_parser(), argparse.ArgumentParser)


def test_paired_flags_preserve_their_destinations():
    arguments = build_parser().parse_args(
        ["electrodynamics", "--dark", "--show", "--name", "example", "results"]
    )

    assert arguments.light is False
    assert arguments.hide is False


def test_help_does_not_import_command_dependencies():
    script = "from md_davis.cli import main; main(['--help']); import sys; assert 'numpy' not in sys.modules"
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, check=False, text=True)

    assert result.returncode == 0, result.stderr


def test_help_lists_supported_commands():
    result = run_cli("--help")

    assert result.returncode == 0
    for command in ("collate", "electrodynamics", "landscape", "sequence", "xvg"):
        assert command in result.stdout


def test_version_uses_package_metadata():
    result = run_cli("--version")

    assert result.returncode == 0
    assert result.stdout.strip() == "md-davis, version 0.4.1"
