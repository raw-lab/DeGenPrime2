#!/usr/bin/env bash
# End-to-end smoke test: unit tests, every CLI subcommand (including the
# built-in MAFFT on unaligned input), and the desktop GUI's HTTP API (session
# token, host check, a full design, the icons). Exits non-zero on any failure,
# so it can run in CI.
set -euo pipefail
cd "$(dirname "$0")/.."
ok()   { printf 'ok   %s\n' "$*"; }
fail() { printf 'FAIL %s\n' "$*" >&2; exit 1; }
BIN=target/release/degenprime
VERSION="$(grep -m1 '^version' Cargo.toml | cut -d'"' -f2)"
TMP="$(mktemp -d)"; GUI_PID=""
trap '[ -n "$GUI_PID" ] && kill "$GUI_PID" 2>/dev/null; rm -rf "$TMP"' EXIT

cargo build --release --locked --quiet && ok "release build"
cargo test --release --locked --quiet >/dev/null 2>&1 && ok "unit tests" || fail "unit tests"
"$BIN" --version | grep -qx "degenprime $VERSION" && ok "version $VERSION" || fail "--version does not say $VERSION"

"$BIN" design -i data/example_anfD.clust -o "$TMP/anfd" --detailed --fasta -q > /dev/null
for f in csv json txt fasta; do [ -s "$TMP/anfd.$f" ] || fail "design did not write .$f"; done
[ "$(wc -l < "$TMP/anfd.csv")" -gt 1 ] || fail "CSV has no primer rows"
ok "design wrote CSV, JSON, report and FASTA"
"$BIN" design -i data/example_anfD.clust -o "$TMP/again" -q > /dev/null
cmp -s "$TMP/anfd.csv" "$TMP/again.csv" && ok "design is deterministic" || fail "two runs differ"
"$BIN" design -t 1 -i data/example_anfD.clust -o "$TMP/one" -q > /dev/null
cmp -s "$TMP/anfd.csv" "$TMP/one.csv" && ok "one thread gives the same result" || fail "--threads 1 differs"
"$BIN" design -i data/example_anfD.clust -o "$TMP/region" --begin 0 --end 700 -q > /dev/null
awk -F, 'NR>1 && ($7 > 700 || $20 > 700) {bad=1} END {exit bad}' "$TMP/region.csv" && ok "design region respected on both strands" || fail "primer outside region"

# Unaligned input: the gaps removed from the bundled MAFFT L-INS-i alignment.
python3 - data/example_anfD.clust "$TMP/genes.fasta" <<'PY'
import sys
seqs = {}
for line in open(sys.argv[1]):
    p = line.split()
    if len(p) == 2 and not line.startswith(("CLUSTAL", " ")):
        seqs[p[0]] = seqs.get(p[0], "") + p[1]
with open(sys.argv[2], "w") as f:
    for name, s in seqs.items():
        f.write(f">{name}\n{s.replace('-', '')}\n")
PY
cp "$TMP/genes.fasta" "$TMP/genes.orig"
"$BIN" align -i "$TMP/genes.fasta" -o "$TMP/realigned.fasta" 2> "$TMP/align.log"
grep -q "L-INS-i" "$TMP/align.log" || fail "align did not choose L-INS-i like mafft --auto"
"$BIN" design -i "$TMP/realigned.fasta" -o "$TMP/realigned" -q > /dev/null
cmp -s "$TMP/anfd.csv" "$TMP/realigned.csv" && ok "align reproduced the MAFFT alignment (same primers)" || fail "built-in MAFFT alignment differs"
"$BIN" design -i "$TMP/genes.fasta" --fasta --save_alignment -q > "$TMP/design.out" 2> "$TMP/design.err"
grep -q "aligned them with MAFFT L-INS-i" "$TMP/design.err" || fail "design did not report the alignment"
cmp -s "$TMP/genes.fasta" "$TMP/genes.orig" || fail "design --fasta overwrote its input"
[ -s "$TMP/genes_primers.fasta" ] && [ -s "$TMP/genes.aligned.fasta" ] || fail "missing genes_primers.fasta or genes.aligned.fasta"
cmp -s "$TMP/genes.csv" "$TMP/anfd.csv" && ok "design aligns unaligned input first, without overwriting it" || fail "unaligned design differs"
if "$BIN" design -i "$TMP/genes.fasta" --align never -o "$TMP/never" -q > /dev/null 2>&1; then :; fi
[ -s "$TMP/never.csv" ] && ok "--align never pads instead" || fail "--align never"

"$BIN" test -p ATGCRYCCGGAATTCCAGTC --partner GACTGGAATTCCGGRYGCAT --json "$TMP/p.json" > "$TMP/test.txt"
grep -q '"fold": 4' "$TMP/p.json" && grep -q "ATGCRYCC" "$TMP/test.txt" && ok "test keeps IUPAC codes (4-fold primer)" || fail "IUPAC handling"
"$BIN" doctor > "$TMP/doctor.txt" || fail "doctor"
grep -q "reproduced the bundled MAFFT alignment" "$TMP/doctor.txt" && ok "doctor (aligner self-test included)" || fail "doctor aligner check"

"$BIN" gui --no_browser --port 0 > "$TMP/gui.log" 2>&1 & GUI_PID=$!
URL=""
for _ in $(seq 1 50); do URL="$(grep -o 'http://127.0.0.1:[0-9]*/' "$TMP/gui.log" | head -1 || true)"; [ -n "$URL" ] && break; sleep 0.1; done
[ -n "$URL" ] || fail "gui did not start"
TOKEN="$(curl -s "$URL" | grep -o '"token":"[0-9a-f]*"' | cut -d'"' -f4)"
[ ${#TOKEN} -eq 32 ] || fail "gui page has no session token"
REQ="$TMP/req.json"; printf '{"input": %s}' "$(python3 -c 'import json,sys;print(json.dumps(open(sys.argv[1]).read()))' "$TMP/genes.fasta")" > "$REQ"
curl -s -X POST -H "X-DGP-Token: $TOKEN" --data-binary @"$REQ" "${URL}api/design" > "$TMP/gui.json"
grep -q '"ok":true' "$TMP/gui.json" && grep -q '"aligned_with":"MAFFT' "$TMP/gui.json" && ok "gui API aligns and designs" || fail "gui API design"
[ "$(curl -s -o /dev/null -w '%{http_code}' -X POST "${URL}api/version")" = 403 ] && ok "gui API rejects calls without the token" || fail "token check"
[ "$(curl -s -o /dev/null -w '%{http_code}' -H 'Host: evil.example' "$URL")" = 403 ] && ok "gui rejects foreign Host headers" || fail "host check"
for p in assets/mark.png assets/logo.webp favicon.ico manifest.webmanifest; do
    [ "$(curl -s -o /dev/null -w '%{http_code}' "$URL$p")" = 200 ] || fail "gui does not serve $p"
done
ok "gui serves the logo, icons and manifest"

for f in web/index.html web/_headers web/manifest.webmanifest web/assets/degenprime.js web/assets/degenprime.css \
         web/assets/degenprime.wasm web/assets/logo.webp web/assets/mark.png web/assets/favicon.ico \
         wordpress/degenprime/degenprime.php assets/icon.ico assets/icon.icns assets/icons/degenprime-256.png \
         assets/windows/degenprime.res assets/windows/degenprime-x86_64-gnu.o src/mafft/LICENSE-MAFFT.txt; do
    [ -s "$f" ] || fail "missing $f"
done
grep -q "FileVersion\", \"$VERSION\"" assets/windows/degenprime.rc || fail "Windows resource is not for $VERSION: run scripts/make_windows_resource.sh"
ok "website, WordPress, icon and licence files present"
printf '\nAll smoke tests passed.\n'
