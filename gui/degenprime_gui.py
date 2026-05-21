#!/usr/bin/env python3
"""
DeGenPrime GUI - Cross-platform PCR primer design with beautiful analytics.

Design goals (informed by the bug review of the previous PyQt6 GUI):
  * Lightweight install: customtkinter + matplotlib + polars + seaborn.
    Works on Linux/macOS/Windows from a single `pip install -r requirements.txt`.
    No PyQt6 / Qt WebEngine binary wheels to fight with.
  * Beautiful outputs: seaborn-rendered embedded plots, polars DataFrame
    backing the results table.
  * Two clear workflows:
      - Design primers from a FASTA/CLUSTAL alignment.
      - Test/analyse an existing primer (your own oligo).
  * Shells out to the Rust `degenprime` binary; reads CSV+JSON back via polars.

Author: Richard Allen White III (rwhit101@uncc.edu)
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox
from typing import Optional

import customtkinter as ctk
import matplotlib
matplotlib.use("Agg")  # noqa: E402  (force backend before pyplot import)
import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# Set seaborn styling once at import time.
sns.set_theme(style="whitegrid", context="talk", palette="viridis")

APP_VERSION = "2.1.0"
DEFAULT_FONT = ("Inter", 12)
HEADER_FONT = ("Inter", 18, "bold")
SUBTITLE_FONT = ("Inter", 13, "bold")


# ---------------------------------------------------------------------------
# Backend binary discovery
# ---------------------------------------------------------------------------

def find_rust_binary() -> Optional[Path]:
    """Locate the compiled `degenprime` binary across platforms.

    Search order:
      1. DEGENPRIME_BIN environment variable
      2. ./target/release/degenprime[.exe] relative to this script
      3. ../target/release/degenprime[.exe] (we live in gui/)
      4. ./bin/degenprime[.exe] (installer drops it here)
      5. PATH (`shutil.which`)
    """
    exe = "degenprime.exe" if platform.system() == "Windows" else "degenprime"
    candidates = []

    env = os.environ.get("DEGENPRIME_BIN")
    if env:
        candidates.append(Path(env))

    script_dir = Path(__file__).resolve().parent
    candidates += [
        script_dir / "target" / "release" / exe,
        script_dir.parent / "target" / "release" / exe,
        script_dir / "bin" / exe,
        script_dir.parent / "bin" / exe,
        Path.cwd() / "target" / "release" / exe,
    ]
    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            return c

    found = shutil.which(exe)
    if found:
        return Path(found)
    return None


# ---------------------------------------------------------------------------
# Worker thread
# ---------------------------------------------------------------------------

class WorkerThread(threading.Thread):
    """Runs the Rust binary off the Tk thread so the UI stays responsive."""

    def __init__(self, binary: Path, args: list[str], on_done):
        super().__init__(daemon=True)
        self.binary = binary
        self.args = args
        self.on_done = on_done

    def run(self):
        try:
            proc = subprocess.run(
                [str(self.binary), *self.args],
                capture_output=True, text=True, timeout=600,
            )
            self.on_done(proc.returncode, proc.stdout, proc.stderr)
        except subprocess.TimeoutExpired:
            self.on_done(-1, "", "Process timed out after 10 minutes")
        except Exception as e:  # noqa: BLE001
            self.on_done(-1, "", f"{type(e).__name__}: {e}")


# ---------------------------------------------------------------------------
# Main application
# ---------------------------------------------------------------------------

class DeGenPrimeApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("DeGenPrime — PCR Primer Design")
        self.geometry("1280x820")
        self.minsize(1000, 700)

        ctk.set_appearance_mode("System")
        ctk.set_default_color_theme("blue")

        self.binary = find_rust_binary()
        self.last_csv: Optional[Path] = None
        self.last_json: Optional[Path] = None
        self.last_df: Optional[pl.DataFrame] = None

        self._build_layout()
        self._check_binary_present()

    # ----- Layout building ------------------------------------------------

    def _build_layout(self):
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self.tabs = ctk.CTkTabview(self, anchor="nw", corner_radius=8)
        self.tabs.grid(row=0, column=0, padx=12, pady=12, sticky="nsew")

        self.tabs.add("Design primers")
        self.tabs.add("Test a primer")
        self.tabs.add("Results")
        self.tabs.add("Plots")
        self.tabs.add("Help")

        self._build_design_tab(self.tabs.tab("Design primers"))
        self._build_test_tab(self.tabs.tab("Test a primer"))
        self._build_results_tab(self.tabs.tab("Results"))
        self._build_plots_tab(self.tabs.tab("Plots"))
        self._build_help_tab(self.tabs.tab("Help"))

        # Status bar
        self.status_var = tk.StringVar(value="Ready.")
        status = ctk.CTkLabel(
            self, textvariable=self.status_var, anchor="w",
            font=DEFAULT_FONT, padx=12,
        )
        status.grid(row=1, column=0, sticky="ew", padx=12, pady=(0, 8))

    def _check_binary_present(self):
        if self.binary is None:
            self.status_var.set(
                "⚠ Rust binary 'degenprime' not found. Build it with: cargo build --release"
            )
        else:
            self.status_var.set(f"Using binary: {self.binary}")

    # ----- Design tab -----------------------------------------------------

    def _build_design_tab(self, parent):
        parent.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(parent, text="Design primers from an alignment",
                             font=HEADER_FONT)
        title.grid(row=0, column=0, sticky="w", padx=16, pady=(16, 4))

        subtitle = ctk.CTkLabel(
            parent,
            text="Drop in a FASTA or CLUSTAL alignment and pick parameters. "
                 "Sensible defaults are loaded — change only what you need.",
            font=DEFAULT_FONT, justify="left", wraplength=900,
        )
        subtitle.grid(row=1, column=0, sticky="w", padx=16, pady=(0, 12))

        # ----- File row -----
        file_frame = ctk.CTkFrame(parent)
        file_frame.grid(row=2, column=0, sticky="ew", padx=16, pady=6)
        file_frame.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(file_frame, text="Input alignment:",
                     font=SUBTITLE_FONT).grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.input_var = tk.StringVar()
        ctk.CTkEntry(file_frame, textvariable=self.input_var,
                     placeholder_text="Select an aligned FASTA / CLUSTAL file…"
                     ).grid(row=0, column=1, sticky="ew", padx=8, pady=10)
        ctk.CTkButton(file_frame, text="Browse…", width=110,
                      command=self._browse_input).grid(row=0, column=2, padx=10, pady=10)

        ctk.CTkLabel(file_frame, text="Output base name:",
                     font=SUBTITLE_FONT).grid(row=1, column=0, padx=10, pady=10, sticky="w")
        self.output_var = tk.StringVar()
        ctk.CTkEntry(file_frame, textvariable=self.output_var,
                     placeholder_text="(optional — auto-derived from input)"
                     ).grid(row=1, column=1, sticky="ew", padx=8, pady=10)
        ctk.CTkButton(file_frame, text="Save as…", width=110,
                      command=self._browse_output).grid(row=1, column=2, padx=10, pady=10)

        # ----- Parameters grid -----
        params_frame = ctk.CTkFrame(parent)
        params_frame.grid(row=3, column=0, sticky="ew", padx=16, pady=12)
        for c in range(4):
            params_frame.grid_columnconfigure(c, weight=1)

        ctk.CTkLabel(params_frame, text="Primer parameters",
                     font=SUBTITLE_FONT).grid(row=0, column=0, columnspan=4,
                                              padx=10, pady=(8, 4), sticky="w")

        # Row 1: primer length + amplicon
        self.min_primer_len = self._spinbox(params_frame, 1, 0, "Min primer length (bp)", 18, 25, 20)
        self.max_primer_len = self._spinbox(params_frame, 1, 1, "Max primer length (bp)", 18, 25, 22)
        self.amplicon = self._spinbox(params_frame, 1, 2, "Min amplicon (bp)", 50, 100000, 100)
        self.max_primers = self._spinbox(params_frame, 1, 3, "Max primer pairs", 1, 50, 5)

        # Row 2: Tm + ΔTm
        self.min_temp = self._spinbox(params_frame, 3, 0, "Min Tm (°C)", 50, 65, 58, step=0.5, float_=True)
        self.max_temp = self._spinbox(params_frame, 3, 1, "Max Tm (°C)", 50, 65, 62, step=0.5, float_=True)
        self.max_tm_diff = self._spinbox(params_frame, 3, 2, "Max ΔTm (°C)", 0, 10, 5, step=0.5, float_=True)
        self.min_gc = self._spinbox(params_frame, 3, 3, "Min GC (%)", 0, 100, 40)

        # Row 3: GC + thermo conditions
        self.max_gc = self._spinbox(params_frame, 5, 0, "Max GC (%)", 0, 100, 60)
        self.primer_conc = self._spinbox(params_frame, 5, 1, "Primer conc (nM)", 50, 5000, 200, step=10, float_=True)
        self.salt_conc = self._spinbox(params_frame, 5, 2, "Na+ conc (mM)", 1, 1000, 50, step=5, float_=True)
        # Empty cell for symmetry

        # Run button
        run_frame = ctk.CTkFrame(parent, fg_color="transparent")
        run_frame.grid(row=4, column=0, sticky="ew", padx=16, pady=(8, 12))
        run_frame.grid_columnconfigure(0, weight=1)

        self.design_progress = ctk.CTkProgressBar(run_frame, mode="indeterminate")
        self.design_progress.grid(row=0, column=0, sticky="ew", padx=8, pady=(0, 8))
        self.design_progress.set(0)

        ctk.CTkButton(
            run_frame, text="Design primers", height=44,
            font=("Inter", 14, "bold"),
            command=self._run_design,
        ).grid(row=1, column=0, sticky="ew", padx=8, pady=4)

    def _spinbox(self, parent, row, col, label, lo, hi, default,
                 step=1, float_=False) -> tk.Variable:
        frame = ctk.CTkFrame(parent, fg_color="transparent")
        frame.grid(row=row, column=col, padx=8, pady=4, sticky="ew")
        frame.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(frame, text=label, font=DEFAULT_FONT, anchor="w").grid(
            row=0, column=0, sticky="ew", padx=4)

        var = tk.DoubleVar(value=float(default)) if float_ else tk.IntVar(value=int(default))
        # CustomTkinter doesn't ship a spinbox; combine entry + +/- buttons.
        ent_frame = ctk.CTkFrame(frame, fg_color="transparent")
        ent_frame.grid(row=1, column=0, sticky="ew")
        ent_frame.grid_columnconfigure(0, weight=1)
        entry = ctk.CTkEntry(ent_frame, textvariable=var, width=80)
        entry.grid(row=0, column=0, sticky="ew", padx=2)

        def adjust(delta):
            try:
                v = (float(var.get()) if float_ else int(var.get())) + delta
                v = max(lo, min(hi, v))
                var.set(v if float_ else int(v))
            except (tk.TclError, ValueError):
                var.set(default)

        ctk.CTkButton(ent_frame, text="−", width=28,
                      command=lambda: adjust(-step)).grid(row=0, column=1, padx=2)
        ctk.CTkButton(ent_frame, text="+", width=28,
                      command=lambda: adjust(step)).grid(row=0, column=2, padx=2)
        return var

    def _browse_input(self):
        p = filedialog.askopenfilename(
            title="Open sequence alignment",
            filetypes=[
                ("FASTA files", "*.fasta *.fa *.faa *.fna *.ffn *.fas"),
                ("CLUSTAL files", "*.clust *.aln *.clustal"),
                ("All files", "*.*"),
            ],
        )
        if p:
            self.input_var.set(p)

    def _browse_output(self):
        p = filedialog.asksaveasfilename(
            title="Save results as",
            defaultextension="",
            filetypes=[("CSV / JSON pair", "*"), ("All files", "*.*")],
        )
        if p:
            # strip extension if user added one — Rust appends .csv/.json
            self.output_var.set(str(Path(p).with_suffix("")))

    # ----- Test tab -------------------------------------------------------

    def _build_test_tab(self, parent):
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(3, weight=1)

        ctk.CTkLabel(parent, text="Test an existing primer",
                     font=HEADER_FONT).grid(row=0, column=0, sticky="w", padx=16, pady=(16, 4))
        ctk.CTkLabel(
            parent,
            text="Paste a primer sequence to compute its Tm, GC%, hairpin score, "
                 "and run it through every quality filter.",
            font=DEFAULT_FONT, justify="left", wraplength=900,
        ).grid(row=1, column=0, sticky="w", padx=16, pady=(0, 12))

        entry_frame = ctk.CTkFrame(parent)
        entry_frame.grid(row=2, column=0, sticky="ew", padx=16, pady=6)
        entry_frame.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(entry_frame, text="Primer sequence (IUPAC accepted):",
                     font=SUBTITLE_FONT).grid(row=0, column=0, sticky="w", padx=10, pady=(8, 0))
        self.test_primer_var = tk.StringVar()
        self.test_entry = ctk.CTkEntry(entry_frame, textvariable=self.test_primer_var,
                                       placeholder_text="e.g. ATGCAATTGGCCGGAATTCC",
                                       font=("JetBrains Mono", 13))
        self.test_entry.grid(row=1, column=0, sticky="ew", padx=10, pady=8)

        # Salt/primer conc inputs
        opt_frame = ctk.CTkFrame(entry_frame, fg_color="transparent")
        opt_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=6)
        for c in range(4):
            opt_frame.grid_columnconfigure(c, weight=1)
        self.test_primer_conc = self._spinbox(opt_frame, 0, 0, "Primer conc (nM)",
                                              50, 5000, 200, step=10, float_=True)
        self.test_salt_conc = self._spinbox(opt_frame, 0, 1, "Na+ conc (mM)",
                                            1, 1000, 50, step=5, float_=True)
        self.test_min_temp = self._spinbox(opt_frame, 0, 2, "Min Tm (°C)", 50, 65, 50, step=0.5, float_=True)
        self.test_max_temp = self._spinbox(opt_frame, 0, 3, "Max Tm (°C)", 50, 65, 65, step=0.5, float_=True)

        ctk.CTkButton(parent, text="Analyse primer", height=40,
                      font=("Inter", 13, "bold"),
                      command=self._run_test).grid(row=3, column=0, sticky="new", padx=16, pady=12)

        # Output area
        self.test_output = ctk.CTkTextbox(parent, font=("JetBrains Mono", 12), height=320)
        self.test_output.grid(row=4, column=0, sticky="nsew", padx=16, pady=(0, 16))
        parent.grid_rowconfigure(4, weight=1)

    # ----- Results tab ----------------------------------------------------

    def _build_results_tab(self, parent):
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(2, weight=1)

        ctk.CTkLabel(parent, text="Designed primer pairs",
                     font=HEADER_FONT).grid(row=0, column=0, sticky="w", padx=16, pady=(16, 4))
        self.results_status = tk.StringVar(
            value="Run a design first — results will appear here."
        )
        ctk.CTkLabel(parent, textvariable=self.results_status,
                     font=DEFAULT_FONT).grid(row=1, column=0, sticky="w", padx=16, pady=(0, 8))

        # Use plain tk.Frame so we can host a Treeview (customtkinter doesn't have a table).
        from tkinter import ttk
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure("Results.Treeview", rowheight=28, font=("Inter", 11),
                        background="#fafafa", fieldbackground="#fafafa")
        style.configure("Results.Treeview.Heading", font=("Inter", 11, "bold"),
                        background="#2b7ce0", foreground="white")

        cols = ("rank", "forward_primer", "fwd_pos", "fwd_tm", "fwd_gc",
                "reverse_primer", "rev_pos", "rev_tm", "rev_gc",
                "amplicon", "tm_diff", "dimer")
        headers = ("#", "Forward", "F pos", "F Tm", "F GC%",
                   "Reverse", "R pos", "R Tm", "R GC%",
                   "Amp (bp)", "ΔTm", "Dimer")
        widths = (40, 200, 60, 60, 60, 200, 60, 60, 60, 80, 60, 60)

        table_wrap = ctk.CTkFrame(parent)
        table_wrap.grid(row=2, column=0, sticky="nsew", padx=16, pady=8)
        table_wrap.grid_rowconfigure(0, weight=1)
        table_wrap.grid_columnconfigure(0, weight=1)

        self.tree = ttk.Treeview(table_wrap, columns=cols, show="headings",
                                 style="Results.Treeview")
        for c, h, w in zip(cols, headers, widths):
            self.tree.heading(c, text=h)
            self.tree.column(c, width=w, anchor="w")
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb = ttk.Scrollbar(table_wrap, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        vsb.grid(row=0, column=1, sticky="ns")

        # Action buttons
        btn_frame = ctk.CTkFrame(parent, fg_color="transparent")
        btn_frame.grid(row=3, column=0, sticky="ew", padx=16, pady=(0, 16))
        ctk.CTkButton(btn_frame, text="Copy selected primer (forward)",
                      command=lambda: self._copy_selected_primer("forward_primer")
                      ).grid(row=0, column=0, padx=4, pady=4)
        ctk.CTkButton(btn_frame, text="Copy selected primer (reverse)",
                      command=lambda: self._copy_selected_primer("reverse_primer")
                      ).grid(row=0, column=1, padx=4, pady=4)
        ctk.CTkButton(btn_frame, text="Open CSV folder",
                      command=self._open_output_folder
                      ).grid(row=0, column=2, padx=4, pady=4)
        ctk.CTkButton(btn_frame, text="Re-export CSV…",
                      command=self._reexport_csv
                      ).grid(row=0, column=3, padx=4, pady=4)

    def _copy_selected_primer(self, col):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("DeGenPrime", "Select a row first.")
            return
        col_idx = {"forward_primer": 1, "reverse_primer": 5}[col]
        value = self.tree.item(sel[0], "values")[col_idx]
        self.clipboard_clear()
        self.clipboard_append(value)
        self.update()
        self.status_var.set(f"Copied {col.replace('_', ' ')}: {value}")

    def _open_output_folder(self):
        if not self.last_csv:
            messagebox.showinfo("DeGenPrime", "No results yet — run a design first.")
            return
        folder = self.last_csv.parent
        try:
            if platform.system() == "Windows":
                os.startfile(folder)  # type: ignore[attr-defined]
            elif platform.system() == "Darwin":
                subprocess.run(["open", str(folder)])
            else:
                subprocess.run(["xdg-open", str(folder)])
        except Exception as e:  # noqa: BLE001
            messagebox.showerror("DeGenPrime", f"Could not open folder: {e}")

    def _reexport_csv(self):
        if self.last_df is None:
            messagebox.showinfo("DeGenPrime", "No results yet — run a design first.")
            return
        p = filedialog.asksaveasfilename(
            title="Export CSV as", defaultextension=".csv",
            filetypes=[("CSV", "*.csv")],
        )
        if not p:
            return
        self.last_df.write_csv(p)
        self.status_var.set(f"Exported {p}")

    # ----- Plots tab ------------------------------------------------------

    def _build_plots_tab(self, parent):
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(1, weight=1)

        ctk.CTkLabel(parent, text="Analytics",
                     font=HEADER_FONT).grid(row=0, column=0, sticky="w", padx=16, pady=(16, 8))

        self.plot_frame = ctk.CTkFrame(parent)
        self.plot_frame.grid(row=1, column=0, sticky="nsew", padx=16, pady=8)
        self.plot_frame.grid_rowconfigure(0, weight=1)
        self.plot_frame.grid_columnconfigure(0, weight=1)

        self.plot_placeholder = ctk.CTkLabel(
            self.plot_frame,
            text="Plots will appear here after a design run.",
            font=DEFAULT_FONT,
        )
        self.plot_placeholder.grid(row=0, column=0)

        self._canvas: Optional[FigureCanvasTkAgg] = None

    # ----- Help tab -------------------------------------------------------

    def _build_help_tab(self, parent):
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        textbox = ctk.CTkTextbox(parent, font=("Inter", 12), wrap="word")
        textbox.grid(row=0, column=0, sticky="nsew", padx=16, pady=16)
        textbox.insert("1.0", HELP_TEXT)
        textbox.configure(state="disabled")

    # ----- Run handlers ---------------------------------------------------

    def _run_design(self):
        if self.binary is None:
            messagebox.showerror("DeGenPrime",
                                 "Rust binary not found.\n\nBuild with: cargo build --release")
            return

        input_path = self.input_var.get().strip()
        if not input_path:
            messagebox.showwarning("DeGenPrime", "Please select an input file.")
            return
        if not Path(input_path).is_file():
            messagebox.showerror("DeGenPrime", f"File not found:\n{input_path}")
            return

        args = ["design",
                "--input", input_path,
                "--min_primer_len", str(int(self.min_primer_len.get())),
                "--max_primer_len", str(int(self.max_primer_len.get())),
                "--amplicon", str(int(self.amplicon.get())),
                "--min_temp", f"{float(self.min_temp.get()):.2f}",
                "--max_temp", f"{float(self.max_temp.get()):.2f}",
                "--min_gc", f"{float(self.min_gc.get()):.2f}",
                "--max_gc", f"{float(self.max_gc.get()):.2f}",
                "--primer_conc", f"{float(self.primer_conc.get()):.2f}",
                "--salt_conc", f"{float(self.salt_conc.get()):.2f}",
                "--max_tm_diff", f"{float(self.max_tm_diff.get()):.2f}",
                "--max_primers", str(int(self.max_primers.get())),
                "--verbose", "--detailed"]
        if self.output_var.get().strip():
            args.extend(["--output", self.output_var.get().strip()])

        self.status_var.set("Designing primers…")
        self.design_progress.start()

        WorkerThread(self.binary, args, self._on_design_done).start()

    def _on_design_done(self, rc, stdout, stderr):
        self.after(0, self._on_design_done_main, rc, stdout, stderr)

    def _on_design_done_main(self, rc, stdout, stderr):
        self.design_progress.stop()
        self.design_progress.set(0)
        if rc != 0:
            msg = stderr.strip() or stdout.strip() or "Unknown error"
            self.status_var.set("Design failed.")
            messagebox.showerror("DeGenPrime — design failed", msg)
            return

        # Parse "CSV: <path>" line from stdout.
        csv_path: Optional[Path] = None
        json_path: Optional[Path] = None
        for line in stdout.splitlines():
            if line.startswith("CSV:"):
                csv_path = Path(line.split(":", 1)[1].strip())
            elif line.startswith("JSON:"):
                json_path = Path(line.split(":", 1)[1].strip())

        if csv_path is None or not csv_path.exists():
            messagebox.showerror("DeGenPrime",
                                 "Design reported success but no CSV was found.")
            return

        self.last_csv = csv_path
        self.last_json = json_path

        # Load with polars.
        df = pl.read_csv(csv_path)
        self.last_df = df
        self._populate_table(df)
        self._draw_plots(df)
        self.results_status.set(
            f"{df.height} primer pair(s) from {self.input_var.get()} "
            f"→ {csv_path.name}"
        )
        self.status_var.set("Design complete.")
        # Auto-switch to results
        self.tabs.set("Results")

    def _populate_table(self, df: pl.DataFrame):
        self.tree.delete(*self.tree.get_children())
        for row in df.iter_rows(named=True):
            self.tree.insert(
                "", "end",
                values=(
                    row["rank"],
                    row["forward_primer"],
                    row["forward_position"],
                    f'{row["forward_tm"]:.2f}',
                    f'{row["forward_gc"]:.1f}',
                    row["reverse_primer"],
                    row["reverse_position"],
                    f'{row["reverse_tm"]:.2f}',
                    f'{row["reverse_gc"]:.1f}',
                    row["amplicon_size"],
                    f'{row["tm_difference"]:.2f}',
                    row["dimer_score"],
                ),
            )

    def _draw_plots(self, df: pl.DataFrame):
        """Render 4 seaborn panels into the Plots tab.

        We feed numpy arrays directly to seaborn rather than converting the
        whole DataFrame to pandas. That keeps polars→pandas (which needs
        pyarrow) off the critical path and is faster for our small tables.
        """
        if self._canvas is not None:
            self._canvas.get_tk_widget().destroy()
        self.plot_placeholder.grid_forget()

        import numpy as np

        rank = df["rank"].to_numpy()
        fwd_tm = df["forward_tm"].to_numpy()
        rev_tm = df["reverse_tm"].to_numpy()
        fwd_gc = df["forward_gc"].to_numpy()
        rev_gc = df["reverse_gc"].to_numpy()
        amp = df["amplicon_size"].to_numpy()
        hp_f = df["hairpin_score_fwd"].to_numpy()
        hp_r = df["hairpin_score_rev"].to_numpy()
        dimer = df["dimer_score"].to_numpy()

        fig = Figure(figsize=(12, 8), dpi=100)
        fig.patch.set_facecolor("white")

        # 1. Forward vs Reverse Tm scatter
        ax1 = fig.add_subplot(2, 2, 1)
        sc = ax1.scatter(fwd_tm, rev_tm, c=amp, s=80 + amp / 20.0,
                         cmap="viridis", edgecolors="white", linewidth=1.5)
        lo = float(min(fwd_tm.min(), rev_tm.min())) - 0.5
        hi = float(max(fwd_tm.max(), rev_tm.max())) + 0.5
        ax1.plot([lo, hi], [lo, hi], "k--", alpha=0.4, label="Tm match")
        ax1.set_title("Forward vs Reverse Tm")
        ax1.set_xlabel("Forward Tm (°C)")
        ax1.set_ylabel("Reverse Tm (°C)")
        cbar = fig.colorbar(sc, ax=ax1)
        cbar.set_label("Amplicon (bp)")
        ax1.legend(loc="best", fontsize=9)

        # 2. Amplicon size by rank
        ax2 = fig.add_subplot(2, 2, 2)
        colors = sns.color_palette("viridis", n_colors=max(len(rank), 1))
        ax2.bar(rank, amp, color=colors)
        ax2.set_title("Amplicon size by rank")
        ax2.set_xlabel("Pair rank")
        ax2.set_ylabel("Amplicon (bp)")
        for r, a in zip(rank, amp):
            ax2.text(int(r), int(a) + max(int(amp.max()), 1) * 0.01,
                     f"{int(a)}", ha="center", va="bottom", fontsize=9)

        # 3. GC content distribution (combined forward + reverse)
        ax3 = fig.add_subplot(2, 2, 3)
        all_gc = np.concatenate([fwd_gc, rev_gc])
        sns.histplot(x=all_gc, bins=8, kde=True, color="#2b7ce0", ax=ax3)
        ax3.set_title("GC content distribution")
        ax3.set_xlabel("GC %")
        ax3.set_ylabel("Count")

        # 4. Secondary-structure heatmap
        ax4 = fig.add_subplot(2, 2, 4)
        matrix = np.array([hp_f, hp_r, dimer], dtype=int)
        sns.heatmap(
            matrix, annot=True, cmap="rocket_r", cbar=False, fmt="d",
            xticklabels=[str(int(r)) for r in rank],
            yticklabels=["Fwd hairpin", "Rev hairpin", "Dimer"],
            ax=ax4,
        )
        ax4.set_title("Secondary-structure scores")
        ax4.set_xlabel("Pair rank")
        ax4.set_ylabel("")

        fig.tight_layout()

        self._canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self._canvas.draw()
        self._canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

    def _run_test(self):
        if self.binary is None:
            messagebox.showerror("DeGenPrime",
                                 "Rust binary not found.\n\nBuild with: cargo build --release")
            return
        primer = self.test_primer_var.get().strip().upper()
        if not primer:
            messagebox.showwarning("DeGenPrime", "Enter a primer sequence first.")
            return

        # Write JSON output to a temp file so we can show the structured result.
        import tempfile
        tmp = Path(tempfile.gettempdir()) / "degenprime_test.json"
        args = ["test",
                "--primer", primer,
                "--json", str(tmp),
                "--primer_conc", f"{float(self.test_primer_conc.get()):.2f}",
                "--salt_conc", f"{float(self.test_salt_conc.get()):.2f}",
                "--min_temp", f"{float(self.test_min_temp.get()):.2f}",
                "--max_temp", f"{float(self.test_max_temp.get()):.2f}"]

        self.status_var.set("Analysing primer…")
        WorkerThread(self.binary, args,
                     lambda rc, out, err: self.after(
                         0, self._on_test_done, rc, out, err, tmp)).start()

    def _on_test_done(self, rc, stdout, stderr, tmp_path):
        self.test_output.delete("1.0", "end")
        if rc != 0:
            self.test_output.insert("1.0", stderr or stdout or "Unknown error.")
            self.status_var.set("Test failed.")
            return

        try:
            data = json.loads(tmp_path.read_text())
        except Exception:  # noqa: BLE001
            self.test_output.insert("1.0", stdout)
            self.status_var.set("Test complete.")
            return

        # Pretty render using polars for filter table.
        passed = all(f["passed"] for f in data.get("filters", []))
        header = "✓ All filters passed" if passed else "✗ Some filters failed"
        text = [
            f"{header}",
            "=" * 60,
            f"Sequence       {data['sequence']}",
            f"Length         {data['length']} bp  ({data['actual_size']} non-gap)",
            f"Tm             {data['tm']:.2f} °C",
            f"GC             {data['gc_content']:.1f} %",
            f"Degenerate     {data['degeneracy_count']} positions",
            f"Gaps           {data['gap_count']}",
            f"Penalty        {data['penalty']:.2f}",
            f"Hairpin score  {data['hairpin_score']}",
            "",
            "Filters:",
        ]
        for f in data.get("filters", []):
            mark = "PASS" if f["passed"] else "FAIL"
            line = f"  [{mark}] {f['name']}"
            if f.get("reason"):
                line += f"  — {f['reason']}"
            text.append(line)

        self.test_output.insert("1.0", "\n".join(text))
        self.status_var.set("Test complete.")


HELP_TEXT = """\
DeGenPrime — PCR primer design from aligned DNA sequences.

WHAT IT DOES
============
Selects the top PCR primer pairs from one or more phylogenetically similar DNA
sequences. Filters: degeneracy, deletions, GC content (incl. 3' end), repeats,
non-complementary ends, melting temperature. Pairs are ranked by ΔTm.

WORKFLOWS
=========
1. Design primers — load an aligned FASTA / CLUSTAL file, tune parameters
   if needed, click "Design primers". Results appear in the Results and
   Plots tabs; CSV and JSON files are written next to the input.

2. Test a primer — paste your own oligo sequence (IUPAC codes allowed),
   click "Analyse primer". Tm, GC, hairpin, and all quality-filter outcomes
   are reported.

INPUT FILES
===========
FASTA: .fasta, .fa, .faa, .fna, .ffn, .fas
CLUSTAL: .clust, .aln, .clustal

Sequences must be aligned (same length). If they aren't, the tool will
pad to the longest with gaps and warn — but for serious work, align with
MAFFT first:

    mafft --localpair --maxiterate 1000 --clustalout input.fa > out.clust

PARAMETER GUIDE
===============
Primer length:    18-25 bp. 20-22 is typical.
Min amplicon:     100 bp by default; raise for longer PCR products.
Tm range:         50-65 °C; default 58-62. Use a tighter range for
                  higher-stringency PCR.
Max ΔTm:          5 °C is generous; 1-2 °C makes a more matched pair.
GC content:       40-60% is the sweet spot for primer stability.
Primer conc:      200 nM for typical reactions.
Salt (Na+):       50 mM matches most PCR buffers; raise for high-salt
                  conditions (e.g. halophile assays).

OUTPUTS
=======
<base>.csv   Tabular primer pairs (rank, sequences, positions, Tm, GC,
             amplicon size, ΔTm, secondary-structure scores).
<base>.json  Full machine-readable report including all settings used.
<base>.txt   Human-friendly summary (only with --detailed).

TIPS
====
* If no primers are found, relax GC or Tm range first, then degeneracy.
* The "Test a primer" tab is great for sanity-checking primers you got
  from another tool or a paper.
* All paths in the CSV are forward-strand positions, including reverse
  primer 5' positions — a fix vs the previous version which reported
  reverse-complement-strand positions.

CITATION
========
Fulghum BW, Tanker S, White III RA (2023). DeGenPrime provides robust
primer design and optimization unlocking the biosphere. bioRxiv.
https://doi.org/10.1101/2023.08.11.553048

CONTACT
=======
Richard Allen White III  <rwhit101@uncc.edu>
https://github.com/raw-lab/DeGenPrime
"""


def main():
    app = DeGenPrimeApp()
    app.mainloop()


if __name__ == "__main__":
    sys.exit(main())
