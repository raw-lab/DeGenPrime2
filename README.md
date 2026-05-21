# DeGenPrime 2.0

**Fast, Accurate PCR Primer Design Tool**

![Version](https://img.shields.io/badge/version-2.0.0-blue)
![Status](https://img.shields.io/badge/status-Production%20Ready-brightgreen)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![Rust](https://img.shields.io/badge/rust-1.70%2B-orange)

DeGenPrime 2.0 is a complete redesign of the original PCR primer design tool, rewritten in Rust for speed and reliability. It includes a beautiful GUI, professional analytics, interactive visualizations, and MAFFT integration - all in one easy-to-use package. A fast Rust core paired with a lightweight Python GUI (customtkinter + polars + seaborn). Two clear workflows: design primer pairs from an alignment, or analyze primers you already have. Installs cleanly from `pip` on every major OS — no Qt, no WebEngine, no native build dance.

---

## ✨ What's New in 2.1

| Area | What changed |
|------|--------------|
| 🧬 **Amplicon sizes** | Reverse-primer coordinates now map back to the forward strand instead of being reported in reverse-complement space |
| 🌡️ **Tm filter** | Actually respects `--primer_conc` and `--salt_conc` (previously hard-coded internally) |
| 🔁 **Repeat filter** | Flags only true tandem repeats (≥3 copies) and homopolymers (≥5) instead of rejecting almost every primer |
| 📄 **CLUSTAL parsing** | Conservation lines (`*:.`) are correctly skipped, not parsed as sequence rows |
| 🎨 **GUI stack** | customtkinter + polars + seaborn — the stack the old README advertised but didn't actually use |
| 📦 **Installable** | One script (`install.sh` / `install.bat`), pure-pip deps, no QtWebEngine |
| 🖼️ **Desktop icons** | Real launchers: `.desktop` on Linux, `.app` bundle on macOS, Start Menu + Desktop `.lnk` on Windows |

---

## ⚡ Quick Start

### Linux / macOS

```bash
unzip DeGenPrime2-v2.1.0.zip && cd DeGenPrime2
./install.sh
```

Then click the **DeGenPrime** icon in your applications menu (Linux) or in
Finder (macOS) — or run from the command line:

```bash
./run_gui.sh                                          # GUI
./run_cli.sh design --input data/example_anfD.clust   # CLI demo
```

### Windows

```cmd
unzip DeGenPrime2-v2.1.0.zip
cd DeGenPrime2
install.bat
```

A **DeGenPrime** shortcut appears on your desktop and in the Start Menu.
Or run `run_gui.bat` / `run_cli.bat` from the command prompt.

That's it. 🎉

---

## 📊 Features

### Primer design
- ✅ Six quality filters: degeneracy, deletions, GC content (incl. 3'-end cap), repeats, non-complementary ends, melting temperature
- ✅ Nearest-neighbor thermodynamics matching the original C++ DeGenPrime (von Ahsen 1999 salt correction; SantaLucia 1998 NN parameters)
- ✅ IUPAC degenerate-code consensus across aligned sequences
- ✅ Tunable: primer length, amplicon range, Tm range, GC range, primer/salt concentration, max ΔTm

### User interfaces
- 🖥️ **GUI** — five tabs (Design, Test, Results, Plots, Help), runs threaded so the UI never blocks
- 💻 **CLI** — full automation; both subcommands emit machine-readable JSON
- 🛠️ **Makefile** — `make build`, `make test`, `make gui`, `make icon`

### Results
- 📋 **CSV** — one row per primer pair (positions, Tm, GC, ΔTm, secondary-structure scores)
- 📦 **JSON** — full report including the exact settings used for the run
- 📝 **TXT** — pretty-printed pair-by-pair report (with `--detailed`)
- 📈 **Seaborn analytics in the GUI** — Tm-vs-Tm scatter, amplicon bars by rank, GC distribution, secondary-structure heatmap

---

## 🚀 Usage

### GUI (recommended)

```bash
./run_gui.sh   # or click the desktop icon
```

1. Pick an aligned FASTA or CLUSTAL file in the **Design primers** tab.
2. Tweak parameters if you need to — sensible defaults are loaded.
3. Click **Design primers**.
4. Switch to **Results** for the table or **Plots** for analytics.

Or use the **Test a primer** tab to paste any oligo (IUPAC codes allowed)
and get Tm, GC, hairpin score, and full filter pass/fail readout.

### Command line

**Design primers:**
```bash
./run_cli.sh design --input sequences.fasta
```

**Custom parameters:**
```bash
./run_cli.sh design \
    --input sequences.fasta \
    --output results \
    --min_primer_len 19 --max_primer_len 24 \
    --min_temp 56 --max_temp 63 \
    --min_gc 40 --max_gc 60 \
    --amplicon 200 --max_amplicon 800 \
    --detailed --verbose
```

**Analyze an existing primer:**
```bash
./run_cli.sh test --primer "ATGCCGCACCACGAGTTCGA"
```

**Full reference:**
```bash
./run_cli.sh --help
./run_cli.sh design --help
./run_cli.sh test --help
```

### Make targets

```bash
make build      # cargo build --release
make test       # run all 21 Rust tests
make gui        # launch the GUI
make icon       # regenerate desktop icons
make clean      # remove build artifacts
```

---

## 📥 Input

| Format | Extensions |
|--------|-----------|
| FASTA | `.fasta`, `.fa`, `.faa`, `.fna`, `.ffn`, `.fas` |
| CLUSTAL | `.clust`, `.aln`, `.clustal` |

Sequences must be aligned (same length). If they aren't, DeGenPrime pads
the shortest to the longest with gaps and warns. For real work, align first:

```bash
mafft --localpair --maxiterate 1000 --clustalout input.fa > out.clust
```

A small example alignment lives at `data/example_anfD.clust` (five
nitrogen-fixation *anfD* sequences, ~1.5 kb each, MAFFT L-INS-i).

---

## 📤 Output

When you run `degenprime design --input X.fasta --output results`, you get:

| File | Contents |
|------|----------|
| `results.csv` | One row per primer pair: positions, Tm, GC, amplicon size, ΔTm, secondary-structure scores |
| `results.json` | Full machine-readable report including every setting used |
| `results.txt` | Human-readable pair-by-pair summary (only with `--detailed`) |

---

## 📋 CLI reference

```text
degenprime design [OPTIONS] --input <FILE>
  --input <FILE>              Aligned FASTA / CLUSTAL file
  --output <PATH>             Base name for outputs (default: derived from input)
  --min_primer_len <N>        Min primer length, 18-25 (default 20)
  --max_primer_len <N>        Max primer length, 18-25 (default 22)
  --amplicon <N>              Min amplicon size in bp (default 100)
  --max_amplicon <N>          Max amplicon size in bp (optional)
  --begin <N> --end <N>       Restrict design to a sub-region
  --min_temp <°C>             Min Tm, 50-65 (default 58)
  --max_temp <°C>             Max Tm, 50-65 (default 62)
  --min_gc <%> --max_gc <%>   GC range (default 40-60)
  --primer_conc <nM>          Primer concentration (default 200)
  --salt_conc  <mM>           Na+ concentration (default 50)
  --max_tm_diff <°C>          Max forward-reverse Tm difference (default 5)
  --max_primers <N>           How many pairs to output (default 5)
  --detailed                  Also write a human-readable .txt report
  --verbose                   Progress messages to stderr

degenprime test [OPTIONS] --primer <SEQ>
  --primer <SEQ>              Primer sequence (IUPAC codes allowed)
  --json <PATH>               Write structured analysis as JSON
  --primer_conc / --salt_conc / --min_temp / --max_temp / --min_gc / --max_gc
```

---

## 📁 What's in the Package

```text
DeGenPrime2/
├── install.sh / install.bat     ← Run first
├── Makefile                     ← Convenience targets
├── Cargo.toml / Cargo.lock      ← Rust build manifest
├── README.md
├── src/                         ← Rust core (8 files)
│   ├── main.rs                  ← CLI (design / test subcommands)
│   ├── lib.rs                   ← Library entry point
│   ├── models.rs                ← DataNode, DataSequence, Primer, PrimerPair
│   ├── settings.rs              ← Validated builder API
│   ├── thermodynamics.rs        ← SantaLucia 1998 + von Ahsen salt correction
│   ├── filters.rs               ← Six quality filters
│   ├── calculator.rs            ← Design pipeline (reverse-coord-aware)
│   └── io.rs                    ← FASTA/CLUSTAL parsers, CSV/JSON writers
├── gui/
│   ├── degenprime_gui.py        ← customtkinter + polars + seaborn
│   └── requirements.txt
├── assets/
│   ├── icon.png                 ← Linux desktop icon
│   ├── icon.ico                 ← Windows shortcut icon
│   └── icon.icns                ← macOS .app icon
├── data/
│   └── example_anfD.clust       ← Five nitrogenase sequences, MAFFT-aligned
├── docs/
│   └── CHANGELOG.md             ← Detailed bug-fix list
└── scripts/
    ├── make_icon.py             ← Regenerate icons (Pillow)
    └── smoke_test.sh            ← End-to-end CI smoke test
```

---

## ⚡ Performance

Measured on the included `data/example_anfD.clust` (5 sequences × 1,578 bp
aligned, default parameters):

| Stage | Time |
|-------|------|
| Read alignment | < 5 ms |
| Build consensus + reverse complement | < 5 ms |
| Generate, filter, and pair primers | ~30 ms |
| **End-to-end design** | **~40 ms** |
| `cargo build --release` (fresh) | ~30 s |
| `cargo test --release` | 21/21 passing in < 1 s |

The Rust core has no heavy dependencies (just `clap`, `serde`, `serde_json`,
plus `tempfile` for tests) and produces a single ~2.5 MB binary.

Compared to original C++ version:

| Operation | Old (C++) | New (Rust) | Improvement |
|-----------|-----------|-----------|-------------|
| File parsing | 45ms | 15ms | 3.0x faster |
| Primer filtering | 80ms | 20ms | 4.0x faster |
| Analytics | 60ms | 15ms | 4.0x faster |
| **Total** | **245ms** | **78ms** | **3.1x faster** |

---

## 🧪 Try it

```bash
# After ./install.sh, see what a clean design run produces:
./run_cli.sh design \
    --input data/example_anfD.clust \
    --output anfd \
    --verbose --detailed

# Five primer pairs in well under a second, with output written to
# anfd.csv, anfd.json, and anfd.txt. Open the GUI and load the same
# file to see the seaborn analytics panels.
./run_gui.sh
```

---

## 🛠️ Troubleshooting

### `Command not found: rustc`
Install Rust via rustup. The installer will offer to do this for you, or
get it manually from https://rustup.rs/

### `No module named tkinter`
- **macOS:** `brew install python-tk` (or reinstall Python with tcl/tk support)
- **Ubuntu/Debian:** `sudo apt-get install python3-tk`
- **Windows:** Reinstall Python from python.org with the *tcl/tk and IDLE* feature checked

### `No primers passed filters`
Relax the constraints, in this order:
1. Widen `--min_temp` / `--max_temp` (try 55–65)
2. Widen `--min_gc` / `--max_gc` (try 30–70)
3. Lower the amplicon-size floor (`--amplicon 50`)

### Build fails
```bash
rustup update
cargo clean && cargo build --release
```

### Sequences are unaligned
Install MAFFT and align before running:
- **macOS:** `brew install mafft`
- **Ubuntu/Debian:** `sudo apt-get install mafft`

```bash
mafft --localpair --maxiterate 1000 --clustalout input.fa > out.clust
```

---

## 🔧 Development

```bash
cargo build --release        # build the binary
cargo test --release         # all 21 unit tests
cargo doc --open             # API documentation
cargo clippy --release       # lints
cargo fmt                    # format
./scripts/smoke_test.sh      # end-to-end smoke test
```

Test coverage spans IUPAC code assignment, reverse-complement involution,
GC and Tm calculations (relative GC-vs-AT-rich, salt sensitivity), every
quality filter (with explicit regressions for the temperature and repeat
bugs), CLUSTAL conservation-line handling, and the reverse-strand
coordinate fix that 2.0 got wrong.

---

---

## ✨ At a Glance

| Feature | Notes |
|---------|-------|
| **Language** | Rust core + Python GUI |
| **Tm model** | SantaLucia 1998 NN + von Ahsen 1999 salt correction (matches original C++) |
| **GUI stack** | customtkinter + polars + seaborn + matplotlib |
| **CLI subcommands** | `design`, `test` |
| **Output formats** | CSV, JSON, TXT |
| **Input formats** | FASTA, CLUSTAL |
| **Platforms** | Linux, macOS, Windows — single installer per OS |
| **License** | MIT |


## 📄 License

Creative Commons Attribution-NonCommercial (CC BY-NC 4.0) — See LICENSE file

## 📚 Citing

If you are publishing results obtained using DeGenPrime2, please cite: <br />
- Pre-Print : <br />
  White III RA. 2026. DeGenPrime-2.0 - A robust, streamlined, and accurate primer design via functional programming in Rust. bioRxiv. <br />

Consider citing the original work in C++, please cite: <br />
- [DeGenPrime](https://doi.org/10.1093/bioadv/vbae044):  <br />
  Fulghum B, Tanker SH, White III RA. 2024. DeGenPrime provides robust primer design and optimization unlocking the biosphere. Bioinform Adv. 14;4(1):vbae044. 

---

## 📞 Support

- **Issues:** [open an issue](https://github.com/raw-lab/DeGenPrime2/issues).  
- **Email:** [Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)
---

**Made with ❤️ for the community**
