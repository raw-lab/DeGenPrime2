//! MAFFT multiple sequence alignment, in pure Rust.
//!
//! This is a vendored copy of **rust-MAFFT 0.2.0** by Lucas Goiriz
//! (<https://github.com/luksgrin/rust-MAFFT>, MIT AND BSD-3-Clause), a port of
//! MAFFT 7.526 by Kazutaka Katoh (BSD-3-Clause). The six crates `mafft-types`,
//! `mafft-scoring`, `mafft-fft`, `mafft-tree`, `mafft-align` and `mafft-core`
//! are included here as modules so that DeGenPrime stays one crate that
//! builds with Rust 1.75 and for WebAssembly. See `LICENSE-MAFFT.txt` in this
//! directory for both licence notices.
//!
//! Changes from the published crates:
//!
//! * crate paths rewritten to module paths (`mafft_core` is
//!   [`pipeline`] here, because a module named `core` would shadow the
//!   standard library's);
//! * edition 2021 instead of 2024 (the code needed no other change) and
//!   dependency versions pinned for Rust 1.75;
//! * `mafft-io` (FASTA files) dropped: DeGenPrime hands sequences over in
//!   memory, so `MafftEngine::align_file` is gone;
//! * floating-point contraction fixed to "never fuse" on every target
//!   (see [`types::fp`]), so an alignment is the same on x86-64, Apple
//!   silicon and in the browser;
//! * memory use bounded for long sequences, marked `degenprime:` in the
//!   code: the profile DP keeps three rows of its score matrix instead of
//!   all of them (bit-identical results), FFT anchors are chained on sparse
//!   rows instead of dense square tables (the same anchors), DP buffers of
//!   unusually large blocks are freed instead of pooled, and a block of
//!   more than [`align::profile::LINEAR_SPACE_CELLS`] cells between anchors
//!   is aligned with the linear-space DP (`MSalignmm`). Twelve 30 kb
//!   genomes now align in about 0.5 GB instead of 2.2 GB;
//! * code examples in doc comments marked as text (they are not doctests of
//!   DeGenPrime's API).
//!
//! Against C MAFFT 7.525 on the bundled anfD genes and on simulated sets of
//! 20-400 sequences of 1-1.6 kb (L-INS-i, FFT-NS-i and FFT-NS-2), this
//! copy gives byte-identical alignments. Two cases can differ, by choosing
//! a different alignment of the same score: a block above
//! `LINEAR_SPACE_CELLS` (C MAFFT runs the full DP up to 30,000 columns; on
//! 50 divergent 10 kb sequences 99.95 % of columns agreed), and sequences
//! longer than 30,000 nt, where C MAFFT switches to `MSalignmm` (more than
//! 99.8 % of columns agreed on 30 kb genomes). DeGenPrime calls this module
//! through [`crate::msa`].

#![allow(
    dead_code,
    unused_imports,
    unused_variables,
    unused_mut,
    unused_assignments,
    unused_parens,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals,
    clippy::all,
    clippy::pedantic,
    rustdoc::all
)]

pub mod align;
pub mod fft;
pub mod pipeline;
pub mod scoring;
pub mod tree;
pub mod types;
