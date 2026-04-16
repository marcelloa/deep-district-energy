@echo off

set "SCRIPT_DIR=%~dp0"

for %%I in ("%SCRIPT_DIR%..\..\bin\DEEP.exe") do set "EXE=%%~fI"
for %%I in ("%SCRIPT_DIR%..\..\examples\NEST_case_study_2\NEST_case_study_2") do set "CASE=%%~fI"
for %%I in ("%SCRIPT_DIR%..\..\data") do set "DATA=%%~fI"

"%EXE%" "%CASE%" "%DATA%" 1
pause