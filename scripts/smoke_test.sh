#!/usr/bin/env bash
# DeGenPrime end-to-end smoke test.
#
# Runs the Rust unit tests, then exercises every CLI surface the GUI
# depends on. Exits non-zero on any failure so it can wire into CI.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

if [ -t 1 ]; then
    GREEN=$'\033[0;32m'; RED=$'\033[0;31m'; BLUE=$'\033[0;34m'
    BOLD=$'\033[1m'; NC=$'\033[0m'
else
    GREEN=''; RED=''; BLUE=''; BOLD=''; NC=''
fi
ok()   { printf "%b✓ %s%b\n"  "$GREEN" "$*" "$NC"; }
say()  { printf "%b▸ %s%b\n"  "$BLUE"  "$*" "$NC"; }
fail() { printf "%b✗ %s%b\n"  "$RED"   "$*" "$NC" >&2; exit 1; }

BIN="$SCRIPT_DIR/../target/release/degenprime"

say "Building release binary"
cargo build --release --quiet
[ -x "$BIN" ] || fail "Binary not built"
ok "Built $BIN"

say "Running Rust unit tests"
cargo test --release --quiet
ok "All Rust tests passed"

# Create a tiny synthetic alignment if no test data is available.
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
cat > "$TMPDIR/sample.fasta" <<'EOF'
>seq1
ATGCGCAATCGCAGCTGTAAGTACCTCGAGGCCATGGAGCTGAAGCGCAACATCGAGAAGCTGAACCGCATCGTGCAGAACATGGAGCGCATCGGCGTGACCTACAACGAGAACGTGATCAAGTACATCGAGTACATCCAGAACAAGCTGGCC
>seq2
ATGCGCAACCGCAGCTGTAAGTACCTGGAGGCCATGGAGCTGAAGCGCAACATCGAGAAGCTGAACCGCATCGTGCAGAACATGGAGCGCATCGGCGTGACCTACAACGAGAACGTGATCAAGTACATCGAGTACATCCAGAACAAGCTGGCC
>seq3
ATGCGCAATCGCAGCTGCAAGTACCTGGAGGCCATGGAGCTGAAGCGCAACATCGAGAAGCTGAACCGCATCGTGCAGAACATGGAGCGCATCGGCGTGACCTACAACGAGAACGTGATCAAGTACATCGAGTACATCCAGAACAAGCTGGCG
EOF

say "Design test (FASTA input)"
"$BIN" design --input "$TMPDIR/sample.fasta" --output "$TMPDIR/result" \
    --amplicon 50 --min_temp 55 --max_temp 65 --min_gc 30 --max_gc 70 \
    --detailed --verbose
[ -f "$TMPDIR/result.csv" ]  || fail "CSV not written"
[ -f "$TMPDIR/result.json" ] || fail "JSON not written"
[ -f "$TMPDIR/result.txt" ]  || fail "TXT not written"
ROWS=$(wc -l < "$TMPDIR/result.csv")
[ "$ROWS" -gt 1 ] || fail "CSV has no primer rows"
ok "design produced $((ROWS-1)) primer pairs"

say "Test a single primer"
"$BIN" test --primer "ATGCCGCACCACGAGTTCGA" --json "$TMPDIR/single.json" \
    --primer_conc 200 --salt_conc 50 >/dev/null
[ -f "$TMPDIR/single.json" ] || fail "Single-primer JSON not written"
ok "test subcommand wrote JSON"

say "Verify JSON contains expected fields"
python3 - "$TMPDIR/result.json" "$TMPDIR/single.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1]))
for k in ("version", "generated_at", "input_file", "template_length",
         "settings", "primer_pairs"):
    assert k in report, f"missing key {k} in design JSON"
assert len(report["primer_pairs"]) > 0
row = report["primer_pairs"][0]
for k in ("rank","forward_primer","reverse_primer","amplicon_size","tm_difference"):
    assert k in row, f"missing key {k} in primer row"

analysis = json.load(open(sys.argv[2]))
for k in ("sequence","length","tm","gc_content","filters"):
    assert k in analysis, f"missing key {k} in analysis JSON"
print("JSON schemas OK")
PY
ok "JSON schemas validated"

# Make sure the GUI script at least parses
say "Verify GUI script syntax"
python3 -c "import ast; ast.parse(open('gui/degenprime_gui.py').read())"
ok "GUI script parses"

printf "\n%b%s%b\n" "$GREEN$BOLD" "All smoke tests passed." "$NC"
