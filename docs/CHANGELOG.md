# Changelog

## 2.3.0 — 2026-10-28

### Renamed
- The project is now just **DeGenPrime**: the repository is `raw-lab/degenprime`, the release zip is
  `degenprime-v2.3.0.zip` with a top folder `degenprime/`, and the program, the installers, the website
  and the WordPress plugin call it DeGenPrime throughout. The command was already `degenprime`.

### Added
- **Built-in MAFFT, in pure Rust.** Unaligned input (sequences of different lengths) is aligned before
  designing, on the command line, in the desktop app and on the website (in the browser), with nothing to
  install. The aligner is a vendored copy of rust-MAFFT 0.2.0 (Lucas Goiriz; MIT and BSD-3-Clause), a port
  of MAFFT 7.526. Both pure-Rust ports we compared gave alignments byte-identical to C MAFFT 7.525 on
  our test sets of up to 400 sequences and 10 kb; rust-MAFFT was chosen because it is a library that
  runs the whole `mafft --auto` pipeline in process, which the WebAssembly website and a single Windows
  executable need, whereas mafft-oxide 0.2.0 replaces MAFFT's helper programs and still needs the
  `mafft` shell script. Strategy `auto` follows `mafft --auto` (L-INS-i, FFT-NS-i or FFT-NS-2 by size);
  FFT-NS-2, FFT-NS-i, L-INS-i, G-INS-i and E-INS-i can be chosen. Sequences on the opposite strand are
  reverse-complemented first, like `mafft --adjustdirection`.
- `degenprime align` writes the alignment; `design` gains `--align auto|always|never`,
  `--align_mode`, `--keep_strand`, `--save_alignment` and `-t/--threads`. The website and desktop app
  have the same options under *Alignment*, and Results offers the alignment for download.
- `degenprime doctor` checks that the aligner reproduces the bundled MAFFT L-INS-i alignment exactly.
- **The DeGenPrime logo** in the app header, the About page and as the icon everywhere: the browser tab,
  the web-app manifest (installable website), the Windows executable (icon and version details, built
  in by `build.rs` from prebuilt resources, so no resource compiler is needed), the Windows Start Menu
  and desktop shortcuts, the macOS app bundle and the Linux menu entry (16 to 512 px). All of them are
  made from `assets/logo-source.png` by `scripts/make_icon.py` and
  `scripts/make_windows_resource.sh` (`make icon`).
- A Cancel button for long designs and alignments (on the website it stops the engine at once; in the
  desktop app the page stops waiting and the computation finishes in the background). Very large
  pasted or opened inputs are kept out of the text box so the page stays responsive.

### Fixed
- **`design --fasta` overwrote its input** when the input was a `.fasta` file and no `-o` was given
  (`genes.fasta` was replaced by the primers). Outputs never overwrite the input now (the primers go to
  `genes_primers.fasta`), and `align` refuses to write over its input.
- **Fewer pairs than asked for.** At most a few thousand best pairs were held; when near-duplicates and
  3′-dimer pairs used them all up, fewer pairs than asked for came back although valid ones existed.
  The search now runs again with a larger pool when that happens.
- **Filter statistics for pairs counted only the pairs the early-exit search happened to visit**, so the
  "why no pairs" hint could blame the wrong rule. Overlap and product-size failures are now counted
  exactly over every forward × reverse combination (in 64-bit counts, which billions of combinations
  need on the 32-bit website), and the hint says when no reverse primer lies downstream of any forward
  primer.
- **CLUSTAL rows with the same name were merged**: names truncated by the aligner can coincide, and two
  sequences became one twice as long. Rows are now joined by their position in each block.
- Files with classic Mac line endings (CR only) were read as one line.
- Settings that are not finite numbers (NaN, infinity from the JSON API or the website) were accepted
  and silently produced nonsense; they are refused with a message. GC limits outside 0–100 %, a
  negative maximum ΔTm and a swapped design region now come with a warning saying what was used.
- Desktop app: one thread per connection without limit (now a pool of 8 workers), uploads were limited
  to 64 MB (now 256 MB), and clients that send `Expect: 100-continue` (curl, scripts) waited a second
  before every upload. The session token is now read from `/dev/urandom` on Linux and macOS instead of
  being derived from hash keys.
- Website: a negative number in a setting reached the engine and came back as a raw parser message (the
  form now refuses negative values and whole numbers below their minimum, naming the field), products
  outside the ladder's range were drawn off the simulated gel, and on phones the manual's parameter
  table made the page wider than the screen (it now shows one card per setting).
- Installers: `install.bat` failed to add DeGenPrime to the PATH of a user who had no user PATH yet,
  could break on user names with quotes, and did not find cargo right after installing Rust;
  `install.sh` stopped silently when it could not ask whether to install Rust. Both build with
  `--locked`, so the pinned dependencies that support Rust 1.75 are always used.

### Changed: speed and memory
- Candidate screening runs in parallel without per-window allocations; Tm and complementarity use
  lookup tables; consensus building is parallel and does not copy the sequences. Results are identical
  to 2.2.0 (checked on a corpus of 21 designs, byte for byte) and do not depend on the thread count.
- Pair search: reverse primers are indexed by position, so each forward primer meets only partners
  that give an allowed product size, in order of cost, and the search stops as soon as no remaining
  partner can beat the pairs held.
- The vendored MAFFT keeps memory bounded on long sequences (three rows of the DP score matrix instead
  of all, sparse FFT-anchor chaining, linear-space DP for blocks over 150 M cells): twelve 30 kb genomes
  align in about 0.5 GB instead of 2.2 GB. The linear-space DP can pick a different alignment of the
  same score than C MAFFT, which keeps the full DP up to 30,000 columns (99.95 % of columns agreed on
  50 divergent 10 kb sequences). See `src/mafft/mod.rs` for every change to the vendored code.

### Housekeeping
- Version 2.3.0; Rust 1.75 is still the minimum (Ubuntu 24.04 apt), checked with `--locked` on 1.75 and
  1.97, for Linux, Windows, macOS and WebAssembly.
- 77 DeGenPrime tests (was 57) plus 165 tests of the vendored aligner; the smoke test covers `align`,
  unaligned input, the input-overwrite fix and the GUI's icons. `cargo fmt` and
  `cargo clippy -D warnings` are clean (`rustfmt.toml`; the vendored aligner keeps its own formatting).
- The crate package lists exactly what it needs to build (`include` in `Cargo.toml`).

## 2.2.0 — 2026-10-07

### Merged
- The two 2.1 archives (`v2_1_0` and `main_best`) had byte-identical Rust
  sources; they differed only in packaging. This release starts from `main` (macOS hardening, LICENSE,
  citations) and brings over `docs/CHANGELOG.md` and `.gitignore`, which only `v2_1_0` had.

### Fixed: engine (each reproduced against 2.1 before fixing, now covered by a test)
- **Nondeterministic results.** Consensus ties were broken by `HashMap` iteration order, which Rust
  randomises per process: six identical runs gave two different primer sets. Ties now resolve in
  A, C, G, T order.
- **IUPAC codes became N** in `degenprime test` (`ATGCRY...` was analysed as `ATGCNN...`, Tm about 8 °C
  low). Codes are kept; ambiguous stacks use the mean of their variants, and the Tm range across
  variants is reported.
- **Design region ignored on the reverse strand.** `--begin/--end` were applied to reverse-complement
  coordinates unchanged, so `--begin 0 --end 300` returned reverse primers at column 1329. The region is
  now mirrored onto the reverse strand.
- **`--begin` alone (or `--end` alone) was silently ignored.**
- **Header-less FASTA was rejected** although format sniffing claimed to accept it; records with an
  empty header were dropped. Both now get generated IDs.
- **Gaps were scored as conserved.** Gaps were excluded before computing conservation, so a column that
  is 80 % gaps scored 100 %, and no column was ever a gap column, so the deletions filter could never
  fire. Gaps now count against conservation and mostly-gap columns are real gap columns.
- **Primer-dimer score checked the wrong ends** (forward 3′ end against the reverse primer's 5′ end):
  a true 3′–3′ dimer scored 0. It now pairs the two 3′ ends.
- **Memory blow-up.** Every forward × reverse combination was materialised before sorting: 2.6 GB for a
  10 kb alignment, killed at 20 kb. A bounded heap with early exit now uses 15 MB at 50 kb.
- **Quadratic consensus building** (`chars().nth(col)` for every cell); now linear.
- **Ranking ignored conservation.** Pairs were sorted by ΔTm with conservation only as a tie-break, so
  every top pair had ΔTm 0.00 regardless of conservation. `Primer`/`PrimerPair` ordering was also
  inconsistent with equality (breaking `Ord`'s contract). See *Changed*.
- Smaller: `sub_seq` could panic on out-of-range input; the Tm code claimed to skip gaps but did not;
  CLUSTAL headers from MUSCLE were parsed as a sequence; doc comments named the wrong salt correction.

### Fixed: packaging
- `Cargo.toml` declared `license = "MIT"` while LICENSE is a non-commercial licence, and the README
  said both. `Cargo.toml` now uses `license-file`; the README states one licence.
- Author e-mail aligned with LICENSE (charlotte.edu).

### Changed
- **One crate for CLI and GUI.** The Python GUI (customtkinter, polars, seaborn) is replaced by
  `degenprime gui`: a standard-library web server on 127.0.0.1 serving the same interface as the website,
  with a per-session token and host-header check. No new dependencies.
- **Ranking** is `score = Pf + Pr + |Tmf − Tmopt| + |Tmr − Tmopt| + |ΔTm|` (conservation first).
  `--rank tm_diff` restores 2.1's ordering.
- **New default filters:** Primer3-style self- and pair complementarity (8 / 3 / 8 / 3) and
  near-duplicate suppression (`--min_3prime_distance 3`).
- CSV/JSON gain `score`, consensus sequences, 1-based `*_from/*_to`, degeneracy fold, conservation and
  complementarity columns; every 2.x column keeps its name and meaning.
- Out-of-range settings are clamped with a warning instead of silently.
- Installers build only the Rust binary and add a launcher for `degenprime gui`.

### Added
- Website (Netlify) and WordPress plugin running the engine as WebAssembly; `scripts/build_wasm.sh`
  builds it with rustup or with Ubuntu's apt Rust 1.75 (std compiled from `rust-src`).
- `degenprime doctor`, `test --partner`, `check` alias, `--mafft`, `--fasta`, `--quiet`, `--rank` and
  options for every new filter; ASCII-art help.
- Filter statistics explaining every run; binding sites of each primer in every sequence.
- 57 unit tests (was 21); `scripts/smoke_test.sh` also exercises the GUI API.

## 2.1.0 — 2026-05-21

### Fixed
- **Build was broken from the start.** `Cargo.toml` declared a nonexistent
  `src/server.rs` binary and pinned an outdated `polars = "0.19"` with a
  feature combination that no longer compiles. Removed both; the Rust core
  now needs only `serde`, `serde_json`, `clap`, and `tempfile`.
- **Amplicon size was nonsense.** `calculator::pair_primers` was computing
  `rev.index - fwd.index + rev.length`, but `rev.index` is a position in
  the reverse-complement coordinate frame. The new `reverse_fwd_strand_position()`
  helper maps it back to the forward strand so amplicon size now matches what
  you'd see on a gel.
- **Temperature filter ignored user settings.** `filters::filter_temperature`
  called a local `calculate_tm` that hard-coded primer concentration (50 nM)
  and salt (50 mM), so the `--primer_conc` and `--salt_conc` flags did nothing.
  The thermodynamics module now takes the runtime settings.
- **Repeat filter killed almost every primer.** The old check
  `w[0]==w[2] && w[1]==w[3]` over a sliding 4-window flagged a single
  isolated `ATAT` substring as a tandem repeat. Now uses correct definitions:
  homopolymer runs ≥5 and tandem dinucleotide repeats with ≥3 copies.
- **CLUSTAL parser swallowed conservation lines.** `io::read_clustal` was
  testing `trimmed.starts_with(' ')` *after* trimming (always false), so the
  conservation track (`*:.` characters) was being parsed as sequence data and
  corrupting the alignment. Fixed to inspect the raw line.
- **`analytics::calculate_temperature_stats` returned the placeholder `60.0`**
  for every primer regardless of sequence. Removed the dead module — real Tm
  comes from `thermodynamics::calculate_tm`.
- **MSA "alignment" was just gap-padding.** `msa::align_pair` right-padded
  the shorter sequence with `-` and claimed to have aligned it. Now if the
  input alignment has ragged lengths we explicitly warn the user and pad to
  the longest sequence (the right answer is to run MAFFT/MUSCLE first; we
  don't pretend to do an MSA in 30 lines of Rust).
- **GUI wasn't what the README claimed.** The shipped GUI was PyQt6 +
  QtWebEngine and never imported `polars` or `seaborn` despite advertising
  them. Rewritten in `customtkinter` (lightweight, pure pip, no Qt build
  dependency), with `polars` loading the CSV directly and `seaborn` rendering
  four embedded plots: Forward vs Reverse Tm, amplicon size by rank, GC
  content distribution, and a secondary-structure score heatmap.
- **GUI passed `--html` to the Rust binary**, which never accepted that flag.
  Removed the broken hand-off; reports are emitted as CSV + JSON + optional
  detailed text.
- **Buttons in the results table were stubs** ("View HTML report", "Export
  CSV") — implemented Copy Forward / Copy Reverse / Open CSV folder /
  Re-export CSV. The CSV is the export now.
- **Install was undocumented and didn't work.** `setup.sh` ran against the
  broken Cargo.toml. New `install.sh` (Linux/macOS) and `install.bat`
  (Windows) check for Rust + Python, build the release binary, set up the
  Python venv, write launcher scripts, and register a desktop icon
  (`.desktop` on Linux, `.app` bundle on macOS, Start Menu + Desktop `.lnk`
  on Windows).
- **No icon existed.** Added `scripts/make_icon.py` (Pillow-based generator)
  and committed `assets/icon.{png,ico,icns}` so all three platforms have
  proper icons without users needing to run anything extra.

### Added
- `degenprime test --primer <SEQ>` subcommand for analyzing a single primer
  (Tm, GC, degeneracy, hairpin, filter pass/fail). The original tool only
  supported the design workflow.
- `--detailed` flag emits a human-readable per-pair report alongside the
  machine-readable CSV/JSON.
- `--verbose` flag prints alignment summary and progress.
- 21 unit tests covering thermodynamics, filters, IO, and the
  reverse-coordinate fix (`cargo test --release`).

### Changed
- Python GUI now uses `customtkinter + polars + seaborn + matplotlib` — the
  stack the README promised. No PyQt, no WebEngine, no electron.
- Settings now have validated builders with sensible defaults (Tm 58–62 °C,
  GC 40–60 %, primer 200 nM, salt 50 mM, dNTPs 0.8 mM).
- Tm formula matches the original C++ DeGenPrime exactly (von Ahsen salt
  correction in the denominator, Ct/4 for non-self-complementary primers),
  so results are comparable to the published tool.

## 2.0.0 — earlier
Initial Rust + PyQt6 attempt. Did not build; GUI was non-functional.
See the bugs list above for what 2.1 inherited and fixed.
