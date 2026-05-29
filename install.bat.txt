@echo off
REM ============================================================
REM DeGenPrime 2.1 - Windows installer
REM ============================================================
REM Builds the Rust backend, sets up a Python venv with GUI deps,
REM and creates Start Menu + Desktop shortcuts (clickable icon).
REM
REM Requires: Rust (https://rustup.rs/) and Python 3.8+ from python.org

setlocal enabledelayedexpansion

set "SCRIPT_DIR=%~dp0"
pushd "%SCRIPT_DIR%"

echo.
echo ============================================================
echo          DeGenPrime 2.1 - Windows installer
echo    Robust PCR primer design with a beautiful Python GUI
echo ============================================================
echo.

REM ---- 1. Rust --------------------------------------------------------
echo [1/5] Checking Rust toolchain...
where cargo >nul 2>&1
if errorlevel 1 (
    echo.
    echo Rust is not installed. Please install Rust first:
    echo     https://rustup.rs/
    echo.
    echo After installing, open a new Command Prompt and re-run install.bat.
    pause
    exit /b 1
)
for /f "tokens=*" %%i in ('rustc --version') do set RUST_VERSION=%%i
echo [OK] %RUST_VERSION%

REM ---- 2. Python ------------------------------------------------------
echo [2/5] Checking Python...
where python >nul 2>&1
if errorlevel 1 (
    echo Python is not installed. Please install Python 3.8+ from:
    echo     https://www.python.org/downloads/
    echo Make sure to check "Add Python to PATH" during install.
    pause
    exit /b 1
)
for /f "tokens=*" %%i in ('python --version') do set PYTHON_VERSION=%%i
echo [OK] %PYTHON_VERSION%

python -c "import tkinter" >nul 2>&1
if errorlevel 1 (
    echo [WARN] tkinter not available in this Python.
    echo        Reinstall Python from python.org with "tcl/tk" feature enabled.
    pause
    exit /b 1
)

REM ---- 3. Build Rust --------------------------------------------------
echo [3/5] Building Rust backend (release mode, may take a few minutes)...
call cargo build --release
if errorlevel 1 (
    echo Failed to build Rust backend.
    pause
    exit /b 1
)
echo [OK] target\release\degenprime.exe built

REM ---- 4. Python venv + deps ------------------------------------------
echo [4/5] Setting up Python environment...
if not exist venv (
    python -m venv venv
)
call venv\Scripts\activate.bat
python -m pip install --quiet --upgrade pip
python -m pip install --quiet -r gui\requirements.txt
if errorlevel 1 (
    echo Failed to install Python dependencies.
    pause
    exit /b 1
)
call venv\Scripts\deactivate.bat
echo [OK] GUI dependencies installed (customtkinter, polars, seaborn, matplotlib)

REM ---- 5. Launchers + Start Menu shortcut ----------------------------
echo [5/5] Creating launchers and Start Menu shortcut...

REM Batch launcher for the GUI
(
    echo @echo off
    echo set "SCRIPT_DIR=%%~dp0"
    echo pushd "%%SCRIPT_DIR%%"
    echo set "DEGENPRIME_BIN=%%SCRIPT_DIR%%target\release\degenprime.exe"
    echo call venv\Scripts\activate.bat
    echo python gui\degenprime_gui.py %%*
    echo popd
) > run_gui.bat

REM Batch launcher for the CLI
(
    echo @echo off
    echo "%%~dp0target\release\degenprime.exe" %%*
) > run_cli.bat

REM VBScript launcher (no console window when double-clicking the shortcut)
(
    echo Set WshShell = CreateObject^("WScript.Shell"^)
    echo WshShell.Run Chr^(34^) ^& "%%~dp0run_gui.bat" ^& Chr^(34^), 0, False
) > run_gui_silent.vbs

echo [OK] Wrote run_gui.bat, run_cli.bat, run_gui_silent.vbs

REM Create Start Menu and Desktop shortcuts using PowerShell.
set "ICON_PATH=%SCRIPT_DIR%assets\icon.ico"
if not exist "%ICON_PATH%" set "ICON_PATH=%SCRIPT_DIR%target\release\degenprime.exe"

powershell -NoProfile -ExecutionPolicy Bypass -Command ^
    "$WS = New-Object -ComObject WScript.Shell;" ^
    "$startMenu = [Environment]::GetFolderPath('Programs');" ^
    "$desktop   = [Environment]::GetFolderPath('Desktop');" ^
    "foreach ($dir in @($startMenu, $desktop)) {" ^
    "  $lnk = $WS.CreateShortcut((Join-Path $dir 'DeGenPrime.lnk'));" ^
    "  $lnk.TargetPath = '%SCRIPT_DIR%run_gui_silent.vbs';" ^
    "  $lnk.WorkingDirectory = '%SCRIPT_DIR%';" ^
    "  $lnk.IconLocation = '%ICON_PATH%,0';" ^
    "  $lnk.Description = 'DeGenPrime - PCR primer design';" ^
    "  $lnk.Save()" ^
    "}"

if errorlevel 1 (
    echo [WARN] Could not create shortcuts automatically.
    echo        You can still launch with: run_gui.bat
) else (
    echo [OK] Start Menu + Desktop shortcuts created.
)

echo.
echo ============================================================
echo                    Installation complete
echo ============================================================
echo.
echo To launch:
echo   GUI:   run_gui.bat   (or click the DeGenPrime shortcut)
echo   CLI:   run_cli.bat design --input your.fasta
echo   Help:  run_cli.bat --help
echo.
popd
pause
