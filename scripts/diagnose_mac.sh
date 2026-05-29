#!/usr/bin/env bash
# DeGenPrime — macOS diagnostic.
#
# Run this when something goes wrong on macOS. It checks every common
# failure point and tells you exactly what to fix.
#
# Usage:  ./scripts/diagnose_mac.sh
#
# Safe to run anywhere; doesn't modify the system.

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ -t 1 ]; then
    GREEN=$'\033[0;32m'; RED=$'\033[0;31m'; YELLOW=$'\033[1;33m'
    BLUE=$'\033[0;34m'; BOLD=$'\033[1m'; NC=$'\033[0m'
else
    GREEN=''; RED=''; YELLOW=''; BLUE=''; BOLD=''; NC=''
fi

pass() { printf "%b✓%b %s\n" "$GREEN" "$NC" "$*"; }
fail() { printf "%b✗%b %s\n" "$RED"   "$NC" "$*"; }
warn() { printf "%b⚠%b %s\n" "$YELLOW" "$NC" "$*"; }
info() { printf "%b▸%b %s\n" "$BLUE"  "$NC" "$*"; }
hdr()  { printf "\n%b== %s ==%b\n" "$BOLD" "$*" "$NC"; }

ISSUES=0
note_issue() { ISSUES=$((ISSUES + 1)); }

# -------------------------------------------------------------------- header
hdr "DeGenPrime macOS Diagnostic"
echo "Project dir: $PROJECT_DIR"
echo "Date:        $(date)"
echo "Shell:       $SHELL"
echo "Uname:       $(uname -a)"

# ------------------------------------------------------------------- OS check
hdr "Operating System"
if [ "$(uname -s)" = "Darwin" ]; then
    pass "Running on macOS $(sw_vers -productVersion 2>/dev/null || echo unknown)"
    pass "Hardware: $(uname -m)"
else
    warn "This script is for macOS. You appear to be on $(uname -s)."
    warn "On Linux, just run ./install.sh and read any errors directly."
    exit 0
fi

# --------------------------------------------------------------- Shell check
hdr "Shell"
info "Login shell:   $SHELL"
info "Current shell: $0 (bash $(bash --version 2>/dev/null | head -1 | awk '{print $4}'))"
if [ -n "${ZSH_VERSION:-}" ]; then
    warn "You appear to have sourced this script into zsh."
    warn "Run it instead with:  ./scripts/diagnose_mac.sh"
fi

# -------------------------------------------------------------- Homebrew
hdr "Homebrew"
if command -v brew >/dev/null 2>&1; then
    pass "brew at $(command -v brew)"
    pass "$(brew --version | head -1)"
else
    warn "Homebrew not installed. We recommend it for managing Python+Tk on macOS."
    info "Install:  /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
fi

# --------------------------------------------------------------- Python
hdr "Python interpreters"
PYTHONS=(
    /opt/homebrew/bin/python3
    /opt/homebrew/bin/python3.13
    /opt/homebrew/bin/python3.12
    /opt/homebrew/bin/python3.11
    /usr/local/bin/python3
    /usr/local/bin/python3.13
    /usr/local/bin/python3.12
    /usr/local/bin/python3.11
    /usr/bin/python3
)
FOUND_WORKING_PY=0
for py in "${PYTHONS[@]}"; do
    [ -x "$py" ] || continue
    VER="$("$py" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])' 2>/dev/null)"
    if [ -z "$VER" ]; then
        warn "$py — present but won't print a version (broken install?)"
        continue
    fi
    # Try to actually create a Tk root (Apple's Tk often imports OK then crashes).
    if "$py" -c 'import tkinter as t; r=t.Tk(); r.destroy()' >/dev/null 2>&1; then
        pass "$py — Python $VER with working Tk"
        FOUND_WORKING_PY=1
    else
        # Try to distinguish "no tkinter at all" from "tkinter imports but Tk fails".
        if "$py" -c 'import tkinter' >/dev/null 2>&1; then
            warn "$py — Python $VER, tkinter imports but Tk() failed"
            warn "    Apple's bundled Tk is deprecated. Use Homebrew Python:"
            warn "    brew install python-tk@${VER%.*}"
        else
            warn "$py — Python $VER, no tkinter module"
            warn "    brew install python-tk@${VER%.*}"
        fi
    fi
done
if [ "$FOUND_WORKING_PY" = 0 ]; then
    fail "No Python with working Tk was found. install.sh will not succeed."
    note_issue
    echo
    echo "    Recommended fix:"
    echo "      brew install python-tk@3.13"
    echo "    Then re-run ./install.sh"
fi

# ---------------------------------------------------------------- Rust
hdr "Rust toolchain"
if command -v cargo >/dev/null 2>&1; then
    pass "$(rustc --version 2>/dev/null || echo "rust present")"
    pass "$(cargo --version)"
else
    warn "Rust not installed. install.sh will offer to install via rustup."
    info "Manual install:  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
fi

# ---------------------------------------------------------- Project state
hdr "Project state"
[ -f "$PROJECT_DIR/Cargo.toml" ] && pass "Cargo.toml present" \
    || { fail "Cargo.toml missing — are you in the right folder?"; note_issue; }

if [ -x "$PROJECT_DIR/target/release/degenprime" ]; then
    pass "Rust binary built: target/release/degenprime"
    VER="$("$PROJECT_DIR/target/release/degenprime" --version 2>/dev/null | head -1)"
    [ -n "$VER" ] && pass "Binary works: $VER" || warn "Binary present but won't print --version"
else
    warn "Rust binary not yet built (run ./install.sh)"
fi

if [ -d "$PROJECT_DIR/venv" ]; then
    pass "Python venv present"
    VENV_PY="$PROJECT_DIR/venv/bin/python"
    if [ -x "$VENV_PY" ]; then
        if "$VENV_PY" -c 'import tkinter as t; r=t.Tk(); r.destroy()' >/dev/null 2>&1; then
            pass "venv Python has working Tk"
        else
            fail "venv Python's Tk is broken — rebuild with: rm -rf venv && ./install.sh"
            note_issue
        fi
        if "$VENV_PY" -c 'import customtkinter, polars, seaborn, matplotlib, pyarrow, numpy' >/dev/null 2>&1; then
            pass "venv has all GUI dependencies"
        else
            fail "venv missing dependencies — rebuild with: ./install.sh"
            note_issue
        fi
    else
        fail "venv/bin/python missing — rebuild with: rm -rf venv && ./install.sh"
        note_issue
    fi
else
    warn "No venv yet (run ./install.sh)"
fi

# --------------------------------------------------------- Launchers
hdr "Launchers"
for f in run_gui.sh run_cli.sh DeGenPrime.app; do
    if [ -e "$PROJECT_DIR/$f" ]; then
        pass "$f present"
    else
        warn "$f missing (re-run ./install.sh)"
    fi
done

# -------------------------------------------------- macOS Gatekeeper / xattr
hdr "Gatekeeper / quarantine"
if [ -d "$PROJECT_DIR/DeGenPrime.app" ]; then
    QUARANTINE="$(xattr -p com.apple.quarantine "$PROJECT_DIR/DeGenPrime.app" 2>/dev/null || true)"
    if [ -n "$QUARANTINE" ]; then
        warn "DeGenPrime.app has quarantine attribute — first launch may be blocked."
        warn "Fix:  xattr -dr com.apple.quarantine '$PROJECT_DIR/DeGenPrime.app'"
    else
        pass "DeGenPrime.app has no quarantine flag"
    fi
fi

# --------------------------------------------------------- Recent log
hdr "Recent GUI log"
LOG="$PROJECT_DIR/gui.log"
if [ -f "$LOG" ]; then
    pass "gui.log present ($(wc -l < "$LOG") lines)"
    info "Last 30 lines:"
    echo "-------- gui.log --------"
    tail -30 "$LOG" | sed 's/^/  /'
    echo "-------------------------"
else
    info "No gui.log yet — try launching the GUI once, then re-run this script."
fi

# ----------------------------------------------------------------- Summary
hdr "Summary"
if [ "$ISSUES" -eq 0 ]; then
    printf "%b%s%b\n" "$GREEN$BOLD" "No critical issues found." "$NC"
    echo
    echo "Try launching:  ./run_gui.sh"
    echo "Or click the DeGenPrime app icon."
else
    printf "%b%s%b\n" "$RED$BOLD" "Found $ISSUES issue(s) — see lines marked ✗ above." "$NC"
    echo
    echo "Most common fix sequence:"
    echo "  1. brew install python-tk@3.13"
    echo "  2. rm -rf venv"
    echo "  3. ./install.sh"
fi
exit "$ISSUES"
