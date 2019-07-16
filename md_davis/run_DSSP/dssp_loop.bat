@echo off
setlocal ENABLEDELAYEDEXPANSION
set OUTPUTDIR=%~dp0
set INPUT=%~1
set SCRIPT_DIR=C:\Users\djmaity\Desktop\DSSP

for /f %%a in ('dir /b %INPUT%') do (
    set filename=%%a
    set PDB=!filename:~0,4!
    echo !PDB!
    "!SCRIPT_DIR!\dssp-2.0.4-win32.exe" -i "%INPUT%\!filename!" -o "!PDB!.dssp"
    )