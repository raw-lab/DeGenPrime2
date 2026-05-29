#!/usr/bin/env bash
# DeGenPrime 2.1 — Linux / macOS installer
#
# What this does:
#   1. Locates a working Python 3.8+ with a *functioning* tkinter
#      (on macOS this often means Homebrew's python-tk@X.Y).
#   2. Locates or installs Rust (with the user's consent).
#   3. Builds the Rust binary in release mode → target/release/degenprime
#   4. Creates a Python venv with the right interpreter and installs GUI deps.
#   5. Drops launcher scripts (run_gui.sh, run_cli.sh) that use absolute
#      paths — so they work when invoked from Finder/Dock with a minimal
#      LaunchServices PATH.
#   6. Installs a desktop icon: .desktop on Linux, .app bundle on macOS.
#
# Explicitly invoked with bash via the shebang, so this runs the same way
# whether your interactive shell is bash, zsh, or fish.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ---- pretty output ---------------------------------------------------------
if [ -t 1 ]; then
    RED=$'\033[0;31m'; GREEN=$'\033[0;32m'; YELLOW=$'\033[1;33m'
    BLUE=$'\033[0;34m'; BOLD=$'\033[1m'; NC=$'\033[0m'
else
    RED=''; GREEN=''; YELLOW=''; BLUE=''; BOLD=''; NC=''
fi
say()   { printf "%b%s%b\n" "$BLUE"   "$*" "$NC"; }
ok()    { printf "%b✓ %s%b\n" "$GREEN"  "$*" "$NC"; }
warn()  { printf "%b⚠ %s%b\n" "$YELLOW" "$*" "$NC"; }
err()   { printf "%b✗ %s%b\n" "$RED"    "$*" "$NC" >&2; }

printf "%b" "$BOLD"
cat <<'BANNER'
+------------------------------------------------------------+
|         DeGenPrime 2.1 — installer (Linux / macOS)         |
|    Robust PCR primer design with a polished Python GUI     |
+------------------------------------------------------------+
BANNER
printf "%b\n" "$NC"

OS="$(uname -s)"
ARCH="$(uname -m)"
case "$OS" in
    Darwin) IS_MAC=true;  say "Detected: macOS ($ARCH)";;
    Linux)  IS_MAC=false; say "Detected: Linux ($ARCH)";;
    *) err "Unsupported OS: $OS (supported: Linux, macOS)"; exit 1;;
esac

# ---- 1. Find a Python with a working tkinter -------------------------------
#
# On macOS, multiple Python installations often coexist:
#   /usr/bin/python3                ← Xcode CLI tools (tk often broken)
#   /opt/homebrew/bin/python3       ← Homebrew on Apple Silicon
#   /usr/local/bin/python3          ← Homebrew on Intel
#   /Library/Frameworks/Python.framework/Versions/X.Y/bin/python3  ← python.org
#
# We search them all and prefer the first one that can actually create a Tk
# window object. That avoids the classic macOS failure where the GUI seems
# to start but immediately crashes because /usr/bin/python3's bundled Tk is
# the deprecated Apple build.

# Candidates to try, in preference order.
PYTHON_CANDIDATES=(
    "${PYTHON3:-}"
    "$(command -v python3 || true)"
    /opt/homebrew/bin/python3
    /opt/homebrew/bin/python3.13
    /opt/homebrew/bin/python3.12
    /opt/homebrew/bin/python3.11
    /usr/local/bin/python3
    /usr/local/bin/python3.13
    /usr/local/bin/python3.12
    /usr/local/bin/python3.11
    /Library/Frameworks/Python.framework/Versions/3.13/bin/python3
    /Library/Frameworks/Python.framework/Versions/3.12/bin/python3
    /Library/Frameworks/Python.framework/Versions/3.11/bin/python3
    /usr/bin/python3
)

check_python() {
    local py="$1"
    [ -n "$py" ] && [ -x "$py" ] || return 1
    # Need 3.8+
    "$py" -c 'import sys; sys.exit(0 if sys.version_info >= (3,8) else 1)' \
        >/dev/null 2>&1 || return 1
    # On macOS we always want the live Tk() test, because Apple's bundled
    # Tk imports cleanly but blows up at Tk() instantiation.
    # On Linux, no DISPLAY = headless install (CI, Docker, server) where
    # Tk() can't be instantiated even though tkinter is fine. In that case
    # we accept import-only as proof.
    if $IS_MAC || [ -n "${DISPLAY:-}" ] || [ -n "${WAYLAND_DISPLAY:-}" ]; then
        "$py" -c 'import tkinter as t; r=t.Tk(); r.destroy()' >/dev/null 2>&1 \
            || return 1
    else
        "$py" -c 'import tkinter' >/dev/null 2>&1 || return 1
    fi
    return 0
}

say "Step 1/5 — finding a Python with a working tkinter"
PYTHON_BIN=""
for py in "${PYTHON_CANDIDATES[@]}"; do
    if check_python "$py"; then
        PYTHON_BIN="$py"
        break
    fi
done

if [ -z "$PYTHON_BIN" ]; then
    err "No Python 3.8+ with a working tkinter was found."
    echo
    echo "macOS tip: install Homebrew Python with tcl-tk support:"
    echo "    brew install python-tk@3.13       # or @3.12 / @3.11"
    echo "    /opt/homebrew/bin/python3 -c 'import tkinter; tkinter.Tk()'"
    echo
    echo "Apple's bundled Tk in /usr/bin/python3 is deprecated and"
    echo "can crash silently on recent macOS versions."
    echo
    echo "Linux tip:"
    echo "    sudo apt-get install -y python3 python3-venv python3-pip python3-tk"
    exit 1
fi
PY_VER="$("$PYTHON_BIN" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')"
ok "Python $PY_VER  ($PYTHON_BIN)"

# ---- 2. Rust ---------------------------------------------------------------
say "Step 2/5 — checking Rust toolchain"
if command -v cargo >/dev/null 2>&1; then
    ok "$(rustc --version 2>/dev/null || echo "rust (cargo $(cargo --version))")"
else
    warn "Rust not found."
    read -r -p "Install Rust now via rustup (recommended)? [y/N] " ans
    if [[ "$ans" =~ ^[Yy]$ ]]; then
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
        # shellcheck disable=SC1091
        source "$HOME/.cargo/env"
        ok "Rust installed."
    else
        err "Cannot continue without Rust. See https://rustup.rs/"
        exit 1
    fi
fi

# ---- 3. Build Rust ---------------------------------------------------------
say "Step 3/5 — building Rust backend (release mode, ~30s)"
cargo build --release
ok "Built target/release/degenprime"

# ---- 4. Python venv + deps -------------------------------------------------
say "Step 4/5 — setting up Python environment"

# If an old venv exists but its python is broken (e.g. user upgraded Homebrew),
# rebuild it from scratch.
VENV_OK=false
if [ -d venv ] && [ -x "venv/bin/python" ]; then
    if $IS_MAC || [ -n "${DISPLAY:-}" ] || [ -n "${WAYLAND_DISPLAY:-}" ]; then
        venv_check='import tkinter as t; r=t.Tk(); r.destroy()'
    else
        venv_check='import tkinter'
    fi
    if venv/bin/python -c "$venv_check" >/dev/null 2>&1; then
        VENV_OK=true
    else
        warn "Existing ./venv has a broken Python — rebuilding"
        rm -rf venv
    fi
fi

if [ "$VENV_OK" = false ]; then
    "$PYTHON_BIN" -m venv venv
    ok "Created ./venv with $PYTHON_BIN"
else
    ok "Existing ./venv reused"
fi

VENV_PY="$SCRIPT_DIR/venv/bin/python"
"$VENV_PY" -m pip install --quiet --upgrade pip
"$VENV_PY" -m pip install --quiet -r gui/requirements.txt
ok "GUI dependencies installed (customtkinter, polars, seaborn, matplotlib, pyarrow, numpy)"

# Sanity-check the venv can import everything before we hand it to the user.
"$VENV_PY" -c 'import customtkinter, polars, seaborn, matplotlib, pyarrow, numpy' \
    >/dev/null 2>&1 || { err "Sanity import failed"; exit 1; }
ok "Sanity import OK"

# ---- 5. Launchers + desktop icon ------------------------------------------
say "Step 5/5 — installing launchers and desktop icon"

# Launchers use ABSOLUTE paths (not just `python` / `degenprime`) so they
# survive launches from Finder/Dock where PATH = /usr/bin:/bin:/usr/sbin:/sbin.
# Errors are tee'd to a logfile in the project dir so users can diagnose
# Finder-launch failures (where stdout/stderr are otherwise invisible).

cat > run_gui.sh <<EOF
#!/usr/bin/env bash
# DeGenPrime GUI launcher (auto-generated by install.sh).
# Safe to invoke from any shell (bash/zsh/fish) and from Finder/Dock.
set -e
DGP_DIR="$SCRIPT_DIR"
DGP_LOG="\$DGP_DIR/gui.log"
cd "\$DGP_DIR"

export DEGENPRIME_BIN="\$DGP_DIR/target/release/degenprime"

# Use the venv python by absolute path so a stripped LaunchServices PATH
# can't break us.
exec "\$DGP_DIR/venv/bin/python" gui/degenprime_gui.py "\$@" \\
    > >(tee -a "\$DGP_LOG") 2> >(tee -a "\$DGP_LOG" >&2)
EOF
chmod +x run_gui.sh
ok "Wrote run_gui.sh"

cat > run_cli.sh <<EOF
#!/usr/bin/env bash
exec "$SCRIPT_DIR/target/release/degenprime" "\$@"
EOF
chmod +x run_cli.sh
ok "Wrote run_cli.sh"

# Desktop icon
if $IS_MAC; then
    APP_DIR="$SCRIPT_DIR/DeGenPrime.app"
    rm -rf "$APP_DIR"
    mkdir -p "$APP_DIR/Contents/MacOS" "$APP_DIR/Contents/Resources"
    cat > "$APP_DIR/Contents/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>DeGenPrime</string>
  <key>CFBundleDisplayName</key><string>DeGenPrime</string>
  <key>CFBundleIdentifier</key><string>edu.uncc.degenprime</string>
  <key>CFBundleVersion</key><string>2.1.0</string>
  <key>CFBundleShortVersionString</key><string>2.1.0</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleSignature</key><string>????</string>
  <key>CFBundleExecutable</key><string>degenprime-launcher</string>
  <key>CFBundleIconFile</key><string>icon.icns</string>
  <key>LSMinimumSystemVersion</key><string>10.14</string>
  <key>NSHighResolutionCapable</key><true/>
  <key>NSRequiresAquaSystemAppearance</key><false/>
</dict>
</plist>
EOF

    # The launcher inside the .app bundle. Critical behaviours:
    #   * absolute path to run_gui.sh so it works regardless of CWD
    #   * any failure shows a native dialog so the user isn't left
    #     wondering why nothing happened
    #   * all output captured to gui.log next to the project
    cat > "$APP_DIR/Contents/MacOS/degenprime-launcher" <<EOF
#!/usr/bin/env bash
DGP_HOME="$SCRIPT_DIR"
LOG="\$DGP_HOME/gui.log"

if [ ! -x "\$DGP_HOME/run_gui.sh" ]; then
    /usr/bin/osascript -e "display dialog \"DeGenPrime cannot find its installation directory at:\\n\\n\$DGP_HOME\\n\\nReinstall by running install.sh in the project folder.\" with title \"DeGenPrime\" buttons {\"OK\"} default button 1 with icon stop"
    exit 1
fi

# Run the GUI and capture its exit status.
"\$DGP_HOME/run_gui.sh" >> "\$LOG" 2>&1 &
disown
EOF
    chmod +x "$APP_DIR/Contents/MacOS/degenprime-launcher"

    if [ -f "$SCRIPT_DIR/assets/icon.icns" ]; then
        cp "$SCRIPT_DIR/assets/icon.icns" "$APP_DIR/Contents/Resources/icon.icns"
    fi

    # Strip any quarantine attribute Gatekeeper may add to a freshly-built
    # bundle (silent if xattr/quarantine isn't present).
    xattr -dr com.apple.quarantine "$APP_DIR" 2>/dev/null || true

    ok "Created DeGenPrime.app"
    echo "    Double-click it, or drag it to /Applications."
    echo "    First-launch tip: if macOS says \"cannot be opened\","
    echo "    right-click the app → Open → Open. That whitelists it."
else
    # Linux .desktop file
    DESKTOP_DIR="$HOME/.local/share/applications"
    ICON_DIR="$HOME/.local/share/icons/hicolor/256x256/apps"
    mkdir -p "$DESKTOP_DIR" "$ICON_DIR"
    DESKTOP_FILE="$DESKTOP_DIR/degenprime.desktop"
    ICON_PATH="$ICON_DIR/degenprime.png"
    if [ -f "$SCRIPT_DIR/assets/icon.png" ]; then
        cp "$SCRIPT_DIR/assets/icon.png" "$ICON_PATH"
    fi
    cat > "$DESKTOP_FILE" <<EOF
[Desktop Entry]
Type=Application
Name=DeGenPrime
GenericName=PCR Primer Designer
Comment=Robust PCR primer design from aligned DNA sequences
Exec=$SCRIPT_DIR/run_gui.sh
Icon=$ICON_PATH
Terminal=false
Categories=Science;Biology;Education;
Keywords=PCR;primer;DNA;bioinformatics;sequence;alignment;
StartupNotify=true
EOF
    chmod +x "$DESKTOP_FILE"
    command -v update-desktop-database >/dev/null 2>&1 && \
        update-desktop-database "$DESKTOP_DIR" >/dev/null 2>&1 || true
    ok "Installed launcher: $DESKTOP_FILE"
    ok "DeGenPrime should appear in your application menu."
fi

# Always end with an explicit success marker so users can see the install
# actually completed (occasionally helpful when the venv pip step is silent).
printf "\n%b%s%b\n\n" "$GREEN$BOLD" "Installation complete." "$NC"
cat <<EOF
To launch:
  GUI:   ${BLUE}./run_gui.sh${NC}    (or click the DeGenPrime icon)
  CLI:   ${BLUE}./run_cli.sh design --input data/example_anfD.clust${NC}
  Help:  ${BLUE}./run_cli.sh --help${NC}

Diagnostics if something goes wrong:
  ${BLUE}./scripts/diagnose_mac.sh${NC}   (macOS only)
  ${BLUE}cat gui.log${NC}                  (last GUI run's stdout/stderr)

Optional: install MAFFT for proper sequence alignment:
$(if $IS_MAC; then echo "  brew install mafft"; else echo "  sudo apt-get install -y mafft"; fi)
EOF
