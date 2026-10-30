#!/usr/bin/env bash
# Build web/assets/degenprime.wasm: the DeGenPrime engine (built-in MAFFT
# included) for the website, the WordPress plugin and the desktop app's web
# interface.
#
#   Rust with the wasm32 standard library installed (rustup users get it with
#     rustup target add wasm32-unknown-unknown, which this script runs)
#   Ubuntu apt Rust (rustc 1.75): sudo apt install rust-src lld-17
#     There is no packaged wasm32 standard library, so this script compiles one
#     from rust-src with -Zbuild-std. No sudo is needed: rust-src is copied to
#     ~/.cache/degenprime together with the matching std lockfile from wasm/.
set -euo pipefail
cd "$(dirname "$0")/.."
TARGET=wasm32-unknown-unknown
OUT=web/assets/degenprime.wasm
ARGS=(rustc --lib --release --locked --target "$TARGET" --no-default-features --crate-type cdylib)
SYSROOT="$(rustc --print sysroot)"

if [ -d "$SYSROOT/lib/rustlib/$TARGET/lib" ]; then
    cargo "${ARGS[@]}"
elif command -v rustup >/dev/null 2>&1; then
    rustup target add "$TARGET"
    cargo "${ARGS[@]}"
else
    SRC="$SYSROOT/lib/rustlib/src/rust"
    [ -d "$SRC/library" ] || { echo "error: Rust source not found. Install it: sudo apt install rust-src" >&2; exit 1; }
    LINKER="${WASM_LD:-}"
    if [ -z "$LINKER" ]; then
        for c in /usr/lib/llvm-17/bin/wasm-ld wasm-ld-17 rust-lld; do
            if command -v "$c" >/dev/null 2>&1; then LINKER="$(command -v "$c")"; break; fi
        done
    fi
    # rustc 1.75's rlibs need an LLVM 17 wasm-ld; LLVM 18+ rejects their metadata.
    [ -n "$LINKER" ] || { echo "error: no wasm linker. Install one: sudo apt install lld-17" >&2; exit 1; }
    RUSTV="$(rustc --version | awk '{print $2}')"
    if [ ! -f "$SRC/Cargo.lock" ] && [ ! -f "$SRC/library/Cargo.lock" ]; then
        CACHE="${XDG_CACHE_HOME:-$HOME/.cache}/degenprime/rust-src-$RUSTV"
        if [ ! -f "$CACHE/Cargo.lock" ]; then
            mkdir -p "$CACHE"
            cp -r "$SRC/library" "$CACHE/"
            if [ -f "wasm/rust-$RUSTV-std.Cargo.lock" ]; then
                cp "wasm/rust-$RUSTV-std.Cargo.lock" "$CACHE/Cargo.lock"
            else
                curl -fsSL "https://raw.githubusercontent.com/rust-lang/rust/$RUSTV/Cargo.lock" -o "$CACHE/Cargo.lock"
            fi
        fi
        export __CARGO_TESTS_ONLY_SRC_ROOT="$CACHE"
    fi
    echo "Building the wasm32 standard library from source with rustc $RUSTV (linker: $LINKER)"
    RUSTC_BOOTSTRAP=1 CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_LINKER="$LINKER" \
        cargo "${ARGS[@]}" -Zbuild-std=std,panic_abort
fi
cp "target/$TARGET/release/degenprime.wasm" "$OUT"
echo "Wrote $OUT ($(wc -c < "$OUT") bytes)"
