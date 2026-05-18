// DeGenPrime - A robust PCR primer design tool using functional Rust patterns
// Rewritten from the original C++ implementation with improved performance and clarity

pub mod models;
pub mod filters;
pub mod calculator;
pub mod io;
pub mod settings;
pub mod thermodynamics;
pub mod analytics;
pub mod visualization;
pub mod msa;

pub use models::{Primer, PrimerPair, DataNode, DataSequence, Sequence};
pub use calculator::PrimerCalculator;
pub use settings::Settings;
pub use analytics::{PrimerAnalytics, AnalyticsEngine};
pub use visualization::VisualizationEngine;
pub use msa::{MSAHandler, MSAMode};

