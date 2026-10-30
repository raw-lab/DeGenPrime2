//! DeGenPrime — degenerate PCR primer design from DNA sequences, aligned or not.
//!
//! One crate, three front ends that share one engine and one JSON API
//! ([`api`]):
//!
//! * the `degenprime` command-line tool (`design`, `align`, `test`, `gui`,
//!   `doctor`);
//! * the desktop app, `degenprime gui`, which serves the web interface on
//!   127.0.0.1 from memory (feature `gui`, standard library only);
//! * a WebAssembly build of this library for the website (Netlify or the
//!   WordPress plugin), compiled with `scripts/build_wasm.sh`.
//!
//! Unaligned input is aligned first with a vendored, pure-Rust port of MAFFT
//! ([`msa`]), so no external aligner is needed on any of them.

pub mod api;
pub mod calculator;
pub mod filters;
#[cfg(all(feature = "gui", not(target_arch = "wasm32")))]
pub mod gui;
pub mod io;
// Vendored rust-MAFFT: keeps its upstream formatting.
#[rustfmt::skip]
mod mafft;
pub mod models;
pub mod msa;
pub mod settings;
pub mod thermodynamics;
#[cfg(target_arch = "wasm32")]
pub mod wasm;

pub use calculator::{DesignOutcome, Explain, PrimerAnalysis, PrimerCalculator, Strand};
pub use filters::FilterResult;
pub use io::{DesignReport, PrimerRow, ResultWriter, SequenceReader};
pub use models::{DataNode, DataSequence, Primer, PrimerPair, Sequence};
pub use settings::{RankBy, Settings};
pub use thermodynamics::{Complementarity, Thermodynamics};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// ASCII-art title used by the CLI help, `doctor` and `gui`.
pub const BANNER: &str = r#"
 ____        ____            ____       _
|  _ \  ___ / ___| ___ _ __ |  _ \ _ __(_)_ __ ___   ___
| | | |/ _ \ |  _ / _ \ '_ \| |_) | '__| | '_ ` _ \ / _ \
| |_| |  __/ |_| |  __/ | | |  __/| |  | | | | | | |  __/
|____/ \___|\____|\___|_| |_|_|   |_|  |_|_| |_| |_|\___|
"#;

/// Banner plus a one-line description with the version.
pub fn banner() -> String {
    format!("{BANNER}  Degenerate PCR primer design, MAFFT built in  |  v{VERSION}\n")
}

/// The anfD example alignment bundled with the desktop build.
#[cfg(not(target_arch = "wasm32"))]
pub const EXAMPLE_ALIGNMENT: &str = include_str!("../data/example_anfD.clust");
