#!/usr/bin/env bash
# DeGenPrime 2.1 — Linux / macOS installer
#
# What this does:
#   1. Checks for Rust + Python 3.8+, installs Rust via rustup if missing
#      (with explicit user consent).
#   2. Builds the Rust binary in release mode → target/release/degenprime
#   3. Sets up a Python venv at ./venv/ and installs GUI deps.
#   4. Creates launch scripts (run_gui.sh, run_cli.sh).
#   5. Optionally installs a desktop icon (Linux .desktop / macOS .app launcher).

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
cat <<'EOF'
+------------------------------------------------------------+
|         DeGenPrime 2.1 — installer (Linux / macOS)         |
|    Robust PCR primer design with a beautiful Python GUI    |
+------------------------------------------------------------+
EOF
printf "%b\n" "$NC"

OS="$(uname -s)"
case "$OS" in
    Darwin) IS_MAC=true;  say "Detected: macOS";;
    Linux)  IS_MAC=false; say "Detected: Linux";;
    *) err "Unsupported OS: $OS (supported: Linux, macOS)"; exit 1;;
esac

# ---- 1. Rust ---------------------------------------------------------------
say "Step 1/5 — checking Rust toolchain"
if command -v cargo >/dev/null 2>&1; then
    ok "$(rustc --version)"
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

# ---- 2. Python -------------------------------------------------------------
say "Step 2/5 — checking Python"
if ! command -v python3 >/dev/null 2>&1; then
    err "python3 not found."
    if $IS_MAC; then echo "  Install with: brew install python@3.12"
    else            echo "  Install with: sudo apt-get install -y python3 python3-venv python3-pip python3-tk"
    fi
    exit 1
fi
PYTHON_VERSION="$(python3 -c 'import sys; print("%d.%d" % sys.version_info[:2])')"
ok "Python $PYTHON_VERSION"

# Make sure tkinter is available — customtkinter needs it.
if ! python3 -c 'import tkinter' >/dev/null 2>&1; then
    warn "tkinter not available."
    if $IS_MAC; then
        echo "  Reinstall Python via Homebrew with tk support: brew install python-tk"
    else
        echo "  Install with: sudo apt-get install -y python3-tk"
    fi
    exit 1
fi

# ---- 3. Build Rust ---------------------------------------------------------
say "Step 3/5 — building Rust backend (release mode)"
cargo build --release
ok "Built target/release/degenprime"

# ---- 4. Python venv + deps -------------------------------------------------
say "Step 4/5 — setting up Python environment"
if [ ! -d "venv" ]; then
    python3 -m venv venv
    ok "Created ./venv"
else
    ok "Existing ./venv reused"
fi

# shellcheck disable=SC1091
source venv/bin/activate
python -m pip install --quiet --upgrade pip
python -m pip install --quiet -r gui/requirements.txt
ok "GUI dependencies installed (customtkinter, polars, seaborn, matplotlib)"
deactivate

# ---- 5. Launchers + desktop icon ------------------------------------------
say "Step 5/5 — installing launchers and desktop icon"

cat > run_gui.sh <<EOF
#!/usr/bin/env bash
cd "$SCRIPT_DIR"
source venv/bin/activate
export DEGENPRIME_BIN="$SCRIPT_DIR/target/release/degenprime"
exec python gui/degenprime_gui.py "\$@"
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
</dict>
</plist>
EOF
    cat > "$APP_DIR/Contents/MacOS/degenprime-launcher" <<EOF
#!/usr/bin/env bash
"$SCRIPT_DIR/run_gui.sh"
EOF
    chmod +x "$APP_DIR/Contents/MacOS/degenprime-launcher"
    # Copy icon if we have one (assets/icon.icns); otherwise the system uses a default.
    if [ -f "$SCRIPT_DIR/assets/icon.icns" ]; then
        cp "$SCRIPT_DIR/assets/icon.icns" "$APP_DIR/Contents/Resources/icon.icns"
    fi
    ok "Created DeGenPrime.app — drag to /Applications or double-click to launch."
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
Exec="$SCRIPT_DIR/run_gui.sh"
Icon=$ICON_PATH
Terminal=false
Categories=Science;Biology;Education;
Keywords=PCR;primer;DNA;bioinformatics;sequence;alignment;
StartupNotify=true
EOF
    chmod +x "$DESKTOP_FILE"
    # Update desktop database if present (silent)
    command -v update-desktop-database >/dev/null 2>&1 && \
        update-desktop-database "$DESKTOP_DIR" >/dev/null 2>&1 || true
    ok "Installed launcher: $DESKTOP_FILE"
    ok "DeGenPrime should now appear in your application menu."
fi

cat <<EOF

${GREEN}${BOLD}Installation complete.${NC}

To launch:
  GUI:   ${BLUE}./run_gui.sh${NC}    (or click the DeGenPrime icon)
  CLI:   ${BLUE}./run_cli.sh design --input your.fasta${NC}
  Help:  ${BLUE}./run_cli.sh --help${NC}

Optional: install MAFFT for proper sequence alignment:
$(if $IS_MAC; then echo "  brew install mafft"; else echo "  sudo apt-get install -y mafft"; fi)
EOF
