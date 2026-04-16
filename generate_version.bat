@echo off
setlocal

set "VERSION="

for /f "delims=" %%i in ('git describe --tags --abbrev=0 2^>nul') do set "VERSION=%%i"

if "%VERSION%"=="" (
    set "VERSION=0.0.0"
)

echo Using version: %VERSION%

set "OUTFILE=%~dp0version.f90"

echo module version_module > "%OUTFILE%"
echo   implicit none >> "%OUTFILE%"
echo   character(len=*), parameter :: app_version = "%VERSION%" >> "%OUTFILE%"
echo end module >> "%OUTFILE%"

echo Generated %OUTFILE%

endlocal