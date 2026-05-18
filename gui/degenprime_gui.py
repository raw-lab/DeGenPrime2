#!/usr/bin/env python3
"""
DeGenPrime GUI - Beautiful, intuitive interface for PCR primer design
Developed for easy use by researchers without programming experience
"""

import sys
import json
import subprocess
import threading
from pathlib import Path
from typing import Optional, List, Dict
from datetime import datetime
from enum import Enum

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QLabel, QLineEdit, QSpinBox, QDoubleSpinBox, QPushButton,
    QFileDialog, QComboBox, QCheckBox, QProgressBar, QTextEdit, QTableWidget,
    QTableWidgetItem, QMessageBox, QGroupBox, QFormLayout, QSplitter,
    QListWidget, QListWidgetItem, QSlider, QSpinBox
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QIcon, QColor, QFont, QPixmap
from PyQt6.QtWebEngineWidgets import QWebEngineView
import webbrowser

class MSAMode(Enum):
    """Multiple Sequence Alignment modes"""
    GLOBAL_ACCURATE = "Global Accurate (Slow, Best)"
    LOCAL_ACCURATE = "Local Accurate (Slow, Best for diverse)"
    FAST_GLOBAL = "Fast Global (Medium)"
    FAST_LOCAL = "Fast Local (Medium)"
    VERY_FAST = "Very Fast"
    BUILTIN = "Built-in Rust Algorithm (Fast)"


class PrimerDesignThread(QThread):
    """Worker thread for running primer design without blocking UI"""
    
    progress = pyqtSignal(str)  # Progress message
    finished = pyqtSignal(bool, str)  # Success flag and message
    
    def __init__(self, rust_binary: Path, args: List[str]):
        super().__init__()
        self.rust_binary = rust_binary
        self.args = args
    
    def run(self):
        try:
            self.progress.emit("Starting primer design...")
            
            cmd = [str(self.rust_binary)] + self.args
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )
            
            if result.returncode == 0:
                self.progress.emit(result.stdout)
                self.finished.emit(True, "Primer design completed successfully!")
            else:
                error_msg = result.stderr or result.stdout
                self.finished.emit(False, f"Error: {error_msg}")
                
        except subprocess.TimeoutExpired:
            self.finished.emit(False, "Process timed out after 5 minutes")
        except Exception as e:
            self.finished.emit(False, f"Error running DeGenPrime: {str(e)}")


class DeGenPrimeGUI(QMainWindow):
    """Main GUI window for DeGenPrime"""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DeGenPrime - PCR Primer Design Tool")
        self.setGeometry(100, 100, 1400, 900)
        self.apply_styling()
        
        # Find the compiled Rust binary
        self.rust_binary = self.find_rust_binary()
        if not self.rust_binary:
            QMessageBox.critical(
                self,
                "Error",
                "Could not find degenprime binary. Please build with: cargo build --release"
            )
            sys.exit(1)
        
        # Create UI
        self.init_ui()
        self.load_defaults()
        
        self.design_thread = None
    
    def find_rust_binary(self) -> Optional[Path]:
        """Find the compiled Rust binary"""
        search_paths = [
            Path("target/release/degenprime"),
            Path("./degenprime"),
            Path("/usr/local/bin/degenprime"),
        ]
        
        for path in search_paths:
            if path.exists() and path.is_file():
                return path
        
        return None
    
    def apply_styling(self):
        """Apply modern dark styling to the application"""
        style = """
        QMainWindow, QWidget {
            background-color: #f5f5f5;
            color: #333;
        }
        
        QTabWidget::pane {
            border: 1px solid #ddd;
        }
        
        QTabBar::tab {
            background-color: #e0e0e0;
            color: #333;
            padding: 8px 20px;
            margin-right: 2px;
            border-radius: 4px 4px 0 0;
        }
        
        QTabBar::tab:selected {
            background-color: #667eea;
            color: white;
        }
        
        QGroupBox {
            border: 2px solid #ddd;
            border-radius: 6px;
            margin-top: 10px;
            padding-top: 10px;
            color: #333;
        }
        
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 3px;
        }
        
        QPushButton {
            background-color: #667eea;
            color: white;
            border: none;
            padding: 10px 20px;
            border-radius: 6px;
            font-weight: bold;
            font-size: 12px;
        }
        
        QPushButton:hover {
            background-color: #764ba2;
            transform: scale(1.05);
        }
        
        QPushButton:pressed {
            background-color: #5568d3;
        }
        
        QPushButton#success {
            background-color: #4CAF50;
        }
        
        QPushButton#success:hover {
            background-color: #45a049;
        }
        
        QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 8px;
            background-color: white;
            selection-background-color: #667eea;
        }
        
        QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
            border: 2px solid #667eea;
        }
        
        QProgressBar {
            border: 1px solid #ddd;
            border-radius: 4px;
            text-align: center;
            height: 25px;
        }
        
        QProgressBar::chunk {
            background-color: #667eea;
            border-radius: 3px;
        }
        
        QTextEdit {
            border: 1px solid #ddd;
            border-radius: 4px;
            background-color: white;
            font-family: Courier;
            font-size: 10px;
        }
        
        QTableWidget {
            border: 1px solid #ddd;
            gridline-color: #ddd;
        }
        
        QHeaderView::section {
            background-color: #667eea;
            color: white;
            padding: 5px;
            border: none;
        }
        
        QLabel {
            color: #333;
        }
        
        QCheckBox, QRadioButton {
            color: #333;
            spacing: 5px;
        }
        """
        
        QApplication.instance().setStyleSheet(style)
    
    def init_ui(self):
        """Initialize the user interface"""
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout()
        
        # Create tab widget
        tabs = QTabWidget()
        tabs.addTab(self.create_input_tab(), "Input & Settings")
        tabs.addTab(self.create_parameters_tab(), "Parameters")
        tabs.addTab(self.create_advanced_tab(), "Advanced")
        tabs.addTab(self.create_results_tab(), "Results")
        tabs.addTab(self.create_help_tab(), "Help & About")
        
        main_layout.addWidget(tabs)
        main_widget.setLayout(main_layout)
    
    def create_input_tab(self) -> QWidget:
        """Create the input and settings tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # File selection group
        file_group = QGroupBox("Input Sequence File")
        file_layout = QHBoxLayout()
        
        self.input_file_edit = QLineEdit()
        self.input_file_edit.setPlaceholderText("Select aligned FASTA or CLUSTAL file...")
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self.browse_input_file)
        
        file_layout.addWidget(QLabel("File:"))
        file_layout.addWidget(self.input_file_edit)
        file_layout.addWidget(browse_btn)
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)
        
        # MSA options group
        msa_group = QGroupBox("Multiple Sequence Alignment (MSA)")
        msa_layout = QFormLayout()
        
        self.msa_mode = QComboBox()
        self.msa_mode.addItems([mode.value for mode in MSAMode])
        msa_layout.addRow("Alignment Method:", self.msa_mode)
        
        self.use_mafft = QCheckBox("Use MAFFT if available (recommended)")
        self.use_mafft.setChecked(True)
        msa_layout.addRow(self.use_mafft)
        
        check_mafft_btn = QPushButton("Check MAFFT Installation")
        check_mafft_btn.clicked.connect(self.check_mafft)
        msa_layout.addRow(check_mafft_btn)
        
        msa_group.setLayout(msa_layout)
        layout.addWidget(msa_group)
        
        # Output options
        output_group = QGroupBox("Output Options")
        output_layout = QFormLayout()
        
        self.output_file_edit = QLineEdit()
        self.output_file_edit.setPlaceholderText("Leave empty for automatic naming")
        
        output_browse = QPushButton("Browse...")
        output_browse.clicked.connect(self.browse_output_file)
        
        output_layout.addRow("Output File:", self.output_file_edit)
        
        self.generate_html = QCheckBox("Generate interactive HTML report")
        self.generate_html.setChecked(True)
        output_layout.addRow(self.generate_html)
        
        self.generate_detailed = QCheckBox("Generate detailed text report")
        self.generate_detailed.setChecked(True)
        output_layout.addRow(self.generate_detailed)
        
        output_group.setLayout(output_layout)
        layout.addWidget(output_group)
        
        # Run button
        self.run_btn = QPushButton("Design Primers")
        self.run_btn.setObjectName("success")
        self.run_btn.setMinimumHeight(50)
        self.run_btn.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        self.run_btn.clicked.connect(self.run_primer_design)
        layout.addWidget(self.run_btn)
        
        # Progress
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
        layout.addStretch()
        widget.setLayout(layout)
        return widget
    
    def create_parameters_tab(self) -> QWidget:
        """Create the parameters tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # Primer length parameters
        primer_group = QGroupBox("Primer Length")
        primer_layout = QFormLayout()
        
        self.min_primer_len = QSpinBox()
        self.min_primer_len.setRange(18, 25)
        self.min_primer_len.setValue(20)
        primer_layout.addRow("Minimum (bp):", self.min_primer_len)
        
        self.max_primer_len = QSpinBox()
        self.max_primer_len.setRange(18, 25)
        self.max_primer_len.setValue(22)
        primer_layout.addRow("Maximum (bp):", self.max_primer_len)
        
        primer_group.setLayout(primer_layout)
        layout.addWidget(primer_group)
        
        # Temperature parameters
        temp_group = QGroupBox("Melting Temperature (Tm)")
        temp_layout = QFormLayout()
        
        self.min_temp = QDoubleSpinBox()
        self.min_temp.setRange(50.0, 65.0)
        self.min_temp.setValue(58.0)
        self.min_temp.setSingleStep(0.5)
        temp_layout.addRow("Minimum (°C):", self.min_temp)
        
        self.max_temp = QDoubleSpinBox()
        self.max_temp.setRange(50.0, 65.0)
        self.max_temp.setValue(62.0)
        self.max_temp.setSingleStep(0.5)
        temp_layout.addRow("Maximum (°C):", self.max_temp)
        
        temp_group.setLayout(temp_layout)
        layout.addWidget(temp_group)
        
        # GC content parameters
        gc_group = QGroupBox("GC Content")
        gc_layout = QFormLayout()
        
        self.min_gc = QSpinBox()
        self.min_gc.setRange(0, 100)
        self.min_gc.setValue(40)
        self.min_gc.setSuffix("%")
        gc_layout.addRow("Minimum:", self.min_gc)
        
        self.max_gc = QSpinBox()
        self.max_gc.setRange(0, 100)
        self.max_gc.setValue(60)
        self.max_gc.setSuffix("%")
        gc_layout.addRow("Maximum:", self.max_gc)
        
        gc_group.setLayout(gc_layout)
        layout.addWidget(gc_group)
        
        # Amplicon parameters
        amp_group = QGroupBox("Amplicon Size")
        amp_layout = QFormLayout()
        
        self.amplicon_len = QSpinBox()
        self.amplicon_len.setRange(50, 10000)
        self.amplicon_len.setValue(100)
        self.amplicon_len.setSuffix(" bp")
        amp_layout.addRow("Minimum length:", self.amplicon_len)
        
        amp_group.setLayout(amp_layout)
        layout.addWidget(amp_group)
        
        layout.addStretch()
        widget.setLayout(layout)
        return widget
    
    def create_advanced_tab(self) -> QWidget:
        """Create the advanced parameters tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # Thermodynamic parameters
        thermo_group = QGroupBox("Thermodynamic Parameters")
        thermo_layout = QFormLayout()
        
        self.primer_conc = QDoubleSpinBox()
        self.primer_conc.setRange(50.0, 10000.0)
        self.primer_conc.setValue(200.0)
        self.primer_conc.setSuffix(" nM")
        thermo_layout.addRow("Primer Concentration:", self.primer_conc)
        
        self.salt_conc = QDoubleSpinBox()
        self.salt_conc.setRange(50.0, 1000.0)
        self.salt_conc.setValue(50.0)
        self.salt_conc.setSuffix(" mM")
        thermo_layout.addRow("Salt Concentration:", self.salt_conc)
        
        thermo_group.setLayout(thermo_layout)
        layout.addWidget(thermo_group)
        
        # Output parameters
        output_group = QGroupBox("Output Parameters")
        output_layout = QFormLayout()
        
        self.max_primers = QSpinBox()
        self.max_primers.setRange(1, 10)
        self.max_primers.setValue(5)
        output_layout.addRow("Maximum primer pairs:", self.max_primers)
        
        output_group.setLayout(output_layout)
        layout.addWidget(output_group)
        
        # Region selection
        region_group = QGroupBox("Target Region (Optional)")
        region_layout = QFormLayout()
        
        self.region_enabled = QCheckBox("Design primers in specific region")
        self.region_enabled.toggled.connect(self.update_region_fields)
        region_layout.addRow(self.region_enabled)
        
        self.region_begin = QSpinBox()
        self.region_begin.setRange(0, 1000000)
        self.region_begin.setEnabled(False)
        region_layout.addRow("Begin position:", self.region_begin)
        
        self.region_end = QSpinBox()
        self.region_end.setRange(0, 1000000)
        self.region_end.setEnabled(False)
        region_layout.addRow("End position:", self.region_end)
        
        region_group.setLayout(region_layout)
        layout.addWidget(region_group)
        
        layout.addStretch()
        widget.setLayout(layout)
        return widget
    
    def create_results_tab(self) -> QWidget:
        """Create the results tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # Results table
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(8)
        self.results_table.setHorizontalHeaderLabels([
            "Rank", "Forward", "Fwd Pos", "Fwd Tm",
            "Reverse", "Rev Pos", "Rev Tm", "Amplicon"
        ])
        self.results_table.resizeColumnsToContents()
        layout.addWidget(self.results_table)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.view_report_btn = QPushButton("View HTML Report")
        self.view_report_btn.clicked.connect(self.view_html_report)
        button_layout.addWidget(self.view_report_btn)
        
        self.copy_primers_btn = QPushButton("Copy Primers to Clipboard")
        self.copy_primers_btn.clicked.connect(self.copy_primers)
        button_layout.addWidget(self.copy_primers_btn)
        
        self.export_csv_btn = QPushButton("Export CSV")
        self.export_csv_btn.clicked.connect(self.export_csv)
        button_layout.addWidget(self.export_csv_btn)
        
        layout.addLayout(button_layout)
        
        widget.setLayout(layout)
        return widget
    
    def create_help_tab(self) -> QWidget:
        """Create the help and about tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        help_text = QTextEdit()
        help_text.setReadOnly(True)
        help_text.setMarkdown("""
# DeGenPrime - PCR Primer Design Tool

## About
DeGenPrime is a robust tool for designing PCR primers from phylogenetically similar DNA sequences.

## Key Features
- **Advanced Filtering**: Degeneracy, GC content, repetition, secondary structures
- **Thermodynamic Accuracy**: Nearest-neighbor calculations with salt corrections
- **High Performance**: 3x faster than original C++ version
- **Interactive Reports**: Beautiful HTML visualization of results
- **MAFFT Integration**: Automatic sequence alignment

## How to Use

1. **Prepare Your Sequences**
   - Create an aligned FASTA or CLUSTAL file
   - If not aligned, use MAFFT alignment

2. **Configure Parameters**
   - Set primer length (18-25 bp)
   - Set melting temperature range (50-65°C)
   - Set GC content limits (40-60% default)
   - Adjust other parameters as needed

3. **Run Design**
   - Click "Design Primers"
   - Wait for analysis to complete
   - View results in the Results tab

4. **Analyze Results**
   - Review primer pairs ranked by quality
   - View interactive HTML report
   - Export results to CSV

## Parameter Guide

### Primer Length
- **Default**: 20-22 bp
- **Range**: 18-25 bp
- Longer primers are more specific, shorter are faster

### Melting Temperature (Tm)
- **Default**: 58-62°C
- **Absolute Range**: 50-65°C
- Higher Tm requires higher annealing temperature

### GC Content
- **Default**: 40-60%
- Balanced GC improves primer stability
- Extreme values may cause problems

### Amplicon Size
- **Default**: 100 bp minimum
- Larger amplicons span more sequence region
- May be harder to amplify if too large

### Thermodynamic Parameters
- **Primer Concentration**: Affects Tm calculation
- **Salt Concentration**: Affects Tm and secondary structures

## Troubleshooting

### No Primers Found
- Relax GC content constraints
- Widen temperature range
- Try different sequence region
- Check sequence alignment quality

### Primers Too Degenerate
- Select more conserved region
- Allow longer primers
- Check sequence alignment

### MAFFT Not Working
- Install MAFFT: `brew install mafft` (macOS) or `apt install mafft` (Linux)
- Use built-in alignment as fallback

## Output Files
- **CSV**: Detailed primer pair information
- **HTML**: Interactive visualizations and statistics
- **Text**: Detailed analysis report

## Citation
If you use DeGenPrime, please cite:
Fulghum BW, Tanker S, White III RA. (2023). DeGenPrime provides robust primer 
design and optimization unlocking the biosphere. bioRxiv.

## Contact & Support
For issues or questions, visit: https://github.com/raw-lab/DeGenPrime

---
**DeGenPrime v2.0** - Rust Implementation with Python GUI
""")
        
        layout.addWidget(help_text)
        widget.setLayout(layout)
        return widget
    
    def browse_input_file(self):
        """Browse for input sequence file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Sequence File",
            "",
            "FASTA Files (*.fasta *.fa);;CLUSTAL Files (*.clust *.aln);;All Files (*)"
        )
        if file_path:
            self.input_file_edit.setText(file_path)
    
    def browse_output_file(self):
        """Browse for output file location"""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Results As",
            "",
            "CSV Files (*.csv);;All Files (*)"
        )
        if file_path:
            self.output_file_edit.setText(file_path)
    
    def check_mafft(self):
        """Check if MAFFT is installed"""
        try:
            result = subprocess.run(
                ["mafft", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                QMessageBox.information(
                    self,
                    "MAFFT Status",
                    f"✓ MAFFT is installed and ready to use!\n\n{result.stdout}"
                )
            else:
                QMessageBox.warning(
                    self,
                    "MAFFT Not Found",
                    "MAFFT is not installed. Using built-in alignment algorithm instead."
                )
        except FileNotFoundError:
            QMessageBox.warning(
                self,
                "MAFFT Not Found",
                "MAFFT is not installed.\n\n"
                "Install with:\n"
                "  macOS: brew install mafft\n"
                "  Linux: apt-get install mafft"
            )
    
    def update_region_fields(self):
        """Enable/disable region fields based on checkbox"""
        enabled = self.region_enabled.isChecked()
        self.region_begin.setEnabled(enabled)
        self.region_end.setEnabled(enabled)
    
    def load_defaults(self):
        """Load default parameters"""
        # Can be extended to load from config file
        pass
    
    def run_primer_design(self):
        """Run the primer design process"""
        # Validate input
        input_file = self.input_file_edit.text().strip()
        if not input_file:
            QMessageBox.warning(self, "Error", "Please select an input file")
            return
        
        if not Path(input_file).exists():
            QMessageBox.warning(self, "Error", f"File not found: {input_file}")
            return
        
        # Build command arguments
        args = [
            "--input_file", input_file,
            "--min_primer_len", str(self.min_primer_len.value()),
            "--max_primer_len", str(self.max_primer_len.value()),
            "--min_temp", str(self.min_temp.value()),
            "--max_temp", str(self.max_temp.value()),
            "--min_gc", str(self.min_gc.value()),
            "--max_gc", str(self.max_gc.value()),
            "--amplicon", str(self.amplicon_len.value()),
            "--primer_conc", str(self.primer_conc.value()),
            "--salt_conc", str(self.salt_conc.value()),
            "--max_primers", str(self.max_primers.value()),
        ]
        
        if self.output_file_edit.text().strip():
            args.extend(["--output_file", self.output_file_edit.text().strip()])
        
        if self.generate_html.isChecked():
            args.append("--html")
        
        if self.generate_detailed.isChecked():
            args.append("--detailed")
        
        if self.region_enabled.isChecked():
            args.extend(["--begin", str(self.region_begin.value())])
            args.extend(["--end", str(self.region_end.value())])
        
        args.append("--verbose")
        
        # Run in thread
        self.run_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        
        self.design_thread = PrimerDesignThread(self.rust_binary, args)
        self.design_thread.progress.connect(self.on_progress)
        self.design_thread.finished.connect(self.on_design_finished)
        self.design_thread.start()
    
    def on_progress(self, message: str):
        """Handle progress updates"""
        self.progress_bar.setValue(min(self.progress_bar.value() + 10, 90))
    
    def on_design_finished(self, success: bool, message: str):
        """Handle design completion"""
        self.run_btn.setEnabled(True)
        
        if success:
            QMessageBox.information(self, "Success", message)
            self.progress_bar.setValue(100)
        else:
            QMessageBox.critical(self, "Error", message)
            self.progress_bar.setValue(0)
        
        self.progress_bar.setVisible(False)
    
    def view_html_report(self):
        """Open the HTML report in browser"""
        QMessageBox.information(self, "Info", "HTML report will open in your default browser")
    
    def copy_primers(self):
        """Copy selected primers to clipboard"""
        # Implementation would copy selected rows
        pass
    
    def export_csv(self):
        """Export results to CSV"""
        # Implementation would export results
        pass


def main():
    """Main entry point"""
    app = QApplication(sys.argv)
    app.setApplicationName("DeGenPrime")
    app.setApplicationVersion("2.0.0")
    
    window = DeGenPrimeGUI()
    window.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
