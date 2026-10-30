#!/usr/bin/env bash
# DeGenPrime installer for Linux and macOS.
#
#   ./install.sh                       build, install to ~/.local/bin, add a launcher
#   PREFIX=/usr/local ./install.sh     install somewhere else (may need sudo)
#
# Needs Rust 1.75 or newer (Ubuntu 24.04: sudo apt install cargo) and nothing
# else: the desktop app is `degenprime gui` and the MAFFT aligner is built in.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
if [ -t 1 ]; then B=$'\033[1m'; G=$'\033[32m'; Y=$'\033[33m'; R=$'\033[31m'; N=$'\033[0m'; else B=''; G=''; Y=''; R=''; N=''; fi
say()  { printf '%s> %s%s\n' "$B" "$*" "$N"; }
ok()   { printf '%sok %s%s\n' "$G" "$*" "$N"; }
warn() { printf '%s!  %s%s\n' "$Y" "$*" "$N"; }
die()  { printf '%sx  %s%s\n' "$R" "$*" "$N" >&2; exit 1; }

PREFIX="${PREFIX:-$HOME/.local}"
BIN_DIR="$PREFIX/bin"
OS="$(uname -s)"
VERSION="$(grep -m1 '^version' Cargo.toml | cut -d'"' -f2)"

say "DeGenPrime $VERSION installer ($OS)"
if [ -x "$HOME/.cargo/bin/cargo" ]; then export PATH="$HOME/.cargo/bin:$PATH"; fi
if ! command -v cargo >/dev/null 2>&1; then
    warn "Rust (cargo) is not installed."
    if command -v apt-get >/dev/null 2>&1; then echo "   Ubuntu/Debian: sudo apt install cargo"; fi
    read -r -p "   Install Rust now with rustup (https://rustup.rs)? [y/N] " ans || ans=""
    case "$ans" in
        [yY]*) curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y; export PATH="$HOME/.cargo/bin:$PATH" ;;
        *) die "Install Rust, then run ./install.sh again." ;;
    esac
fi
RUSTV="$(rustc --version | awk '{print $2}')"
if [ "$(printf '%s\n1.75.0\n' "$RUSTV" | sort -V | head -1)" != "1.75.0" ]; then
    die "rustc $RUSTV is too old; DeGenPrime needs 1.75 or newer."
fi
ok "rustc $RUSTV"

say "Building (release)"
cargo build --release --locked
mkdir -p "$BIN_DIR"
install -m 755 target/release/degenprime "$BIN_DIR/degenprime"
ok "Installed $BIN_DIR/degenprime"
case ":$PATH:" in
    *":$BIN_DIR:"*) ;;
    *) warn "$BIN_DIR is not on your PATH. Add to your shell profile:  export PATH=\"$BIN_DIR:\$PATH\"" ;;
esac

# Desktop launchers open a terminal window running `degenprime gui`, so the
# app is visible while it runs and stops when that window is closed.
if [ "$OS" = "Linux" ]; then
    DATA="${XDG_DATA_HOME:-$HOME/.local/share}"
    mkdir -p "$DATA/applications"
    for s in 16 24 32 48 64 128 256 512; do
        mkdir -p "$DATA/icons/hicolor/${s}x${s}/apps"
        cp "assets/icons/degenprime-$s.png" "$DATA/icons/hicolor/${s}x${s}/apps/degenprime.png"
    done
    touch "$DATA/icons/hicolor"  # tells the desktop to look for new icons
    cat > "$DATA/applications/degenprime.desktop" <<EOF
[Desktop Entry]
Type=Application
Version=1.0
Name=DeGenPrime
GenericName=Primer design
Comment=Degenerate PCR primer design
Keywords=PCR;primer;degenerate;DNA;alignment;MAFFT;bioinformatics;
TryExec=$BIN_DIR/degenprime
Exec=$BIN_DIR/degenprime gui
Icon=degenprime
Terminal=true
Categories=Science;Biology;
EOF
    if command -v update-desktop-database >/dev/null 2>&1; then update-desktop-database -q "$DATA/applications" 2>/dev/null || true; fi
    ok "Added DeGenPrime to your applications menu"
elif [ "$OS" = "Darwin" ]; then
    APP="$HOME/Applications/DeGenPrime.app"
    mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
    cp assets/icon.icns "$APP/Contents/Resources/DeGenPrime.icns"
    cat > "$APP/Contents/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>CFBundleName</key><string>DeGenPrime</string>
  <key>CFBundleDisplayName</key><string>DeGenPrime</string>
  <key>CFBundleIdentifier</key><string>org.rawlab.degenprime</string>
  <key>CFBundleVersion</key><string>$VERSION</string>
  <key>CFBundleShortVersionString</key><string>$VERSION</string>
  <key>CFBundleExecutable</key><string>DeGenPrime</string>
  <key>CFBundleIconFile</key><string>DeGenPrime</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>NSHighResolutionCapable</key><true/>
</dict></plist>
EOF
    # Absolute paths: Finder and the Dock start apps with a minimal PATH.
    printf '#!/bin/bash\nexec "%s/degenprime" gui\n' "$BIN_DIR" > "$APP/Contents/Resources/degenprime-gui.command"
    printf '#!/bin/bash\nexec /usr/bin/open -a Terminal "%s/Contents/Resources/degenprime-gui.command"\n' "$APP" > "$APP/Contents/MacOS/DeGenPrime"
    chmod +x "$APP/Contents/Resources/degenprime-gui.command" "$APP/Contents/MacOS/DeGenPrime"
    xattr -dr com.apple.quarantine "$APP" 2>/dev/null || true
    touch "$APP"  # makes Finder pick up a changed icon
    ok "Created $APP"
fi

"$BIN_DIR/degenprime" doctor || warn "doctor reported a problem (see above)"
echo
ok "Done. Desktop app: degenprime gui    Command line: degenprime --help"
