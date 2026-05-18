@echo off
REM DeGenPrime 2.0 - Setup Script for Windows

setlocal enabledelayedexpansion

echo.
echo ============================================================
echo          DeGenPrime 2.0 - Setup Script
echo    Fast, Accurate PCR Primer Design Tool
echo ============================================================
echo.

REM Check if Rust is installed
echo [1/5] Checking Rust installation...
rustc --version >nul 2>&1
if !errorlevel! neq 0 (
    echo.
    echo Rust is not installed. Please install from:
    echo https://rustup.rs/
    echo.
    pause
    exit /b 1
) else (
    for /f "tokens=*" %%i in ('rustc --version') do set RUST_VERSION=%%i
    echo [OK] !RUST_VERSION!
)

REM Check if Python is installed
echo [2/5] Checking Python installation...
python --version >nul 2>&1
if !errorlevel! neq 0 (
    echo.
    echo Python is not installed. Please install from:
    echo https://www.python.org/downloads/
    echo.
    pause
    exit /b 1
) else (
    for /f "tokens=*" %%i in ('python --version') do set PYTHON_VERSION=%%i
    echo [OK] !PYTHON_VERSION!
)

REM Check for MAFFT
echo [3/5] Checking for MAFFT (optional)...
mafft --version >nul 2>&1
if !errorlevel! neq 0 (
    echo [SKIP] MAFFT not found - built-in fallback available
    echo        To install MAFFT, download from:
    echo        https://mafft.cbrc.jp/alignment/software/
) else (
    for /f "tokens=*" %%i in ('mafft --version ^| findstr /i mafft') do set MAFFT_VERSION=%%i
    echo [OK] MAFFT found
)

REM Build Rust backend
echo [4/5] Building Rust backend...
echo This may take 2-3 minutes on first build...
call cargo build --release
if !errorlevel! neq 0 (
    echo.
    echo Failed to build Rust backend
    pause
    exit /b 1
)
echo [OK] Rust backend built

REM Install Python dependencies
echo [5/5] Installing Python GUI dependencies...
pip install --upgrade pip
pip install -r gui\requirements.txt
if !errorlevel! neq 0 (
    echo.
    echo Failed to install Python dependencies
    pause
    exit /b 1
)
echo [OK] Python dependencies installed

REM Create batch scripts
echo Creating launch scripts...

(
    echo @echo off
    echo python gui\degenprime_gui.py %%*
) > run_gui.bat

(
    echo @echo off
    echo target\release\degenprime.exe %%*
) > run_cli.bat

REM Final message
echo.
echo ============================================================
echo             Setup Complete!
echo ============================================================
echo.
echo Quick start:
echo.
echo   GUI (recommended):
echo     run_gui.bat
echo.
echo   Command Line:
echo     run_cli.bat --input_file sequences.fasta
echo.
echo   Get help:
echo     target\release\degenprime.exe --help
echo.
echo Documentation:
echo   - INSTALL.md: Setup instructions
echo   - EXAMPLES.md: Usage examples
echo   - COMPLETE_SUMMARY.md: Feature overview
echo.
pause
