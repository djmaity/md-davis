#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from cx_Freeze import setup, Executable


# Dependencies are automatically detected, but it might need fine tuning.
# "packages": ["os"] is used as example only
build_exe_options = {"packages": ["os"], "excludes": ["tkinter"]}

# base="Win32GUI" should be used only for Windows GUI app
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(
    name = "md-davis",
    version = "0.3.0",
    description='A tool for comparative analysis of molecular dynamics '
                'simulations of proteins.',
    options = {"build_exe": build_exe_options},
    executables = [Executable("md-davis.py", base=base)]
)
