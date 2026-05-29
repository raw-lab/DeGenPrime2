//! DeGenPrime 2.1 — robust PCR primer design for phylogenetically similar DNA.
//!
//! Rust core library. The binary `degenprime` is a thin CLI on top of this
//! API; the Python GUI shells out to that binary and reads the CSV / JSON
//! output back in.

pub mod calculator;
pub mod filters;
pub mod io;
pub mod models;
pub mod settings;
pub mod thermodynamics;

pub use calculator::{PrimerAnalysis, PrimerCalculator};
pub use filters::FilterResult;
pub use io::{DesignReport, PrimerRow, ResultWriter, SequenceReader};
pub use models::{DataNode, DataSequence, Primer, PrimerPair, Sequence};
pub use settings::Settings;
pub use thermodynamics::Thermodynamics;
