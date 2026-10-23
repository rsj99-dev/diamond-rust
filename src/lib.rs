//! A Rust reimplementation of [DIAMOND](https://github.com/bbuchfink/diamond).
//!
//! The pipeline mirrors the C++ program:
//!
//! 1. [`db`] builds a compact binary database; the seed index is created in
//!    memory at search time, exactly like DIAMOND 2.x.
//! 2. [`seed`] hashes reduced-alphabet spaced seeds from both the query and
//!    the reference and joins them (`index_chunks` limits peak memory).
//! 3. [`ungapped`] extends every hit with the seed-diagonal X-drop algorithm;
//!    sub-optimal hits on the same diagonal are dropped.
//! 4. [`gapped`] runs a banded affine-gap Smith-Waterman with traceback.
//! 5. [`search`] turns raw scores into bit scores/E-values, sorts and culls.
//! 6. [`output`] prints BLAST tabular (outfmt 6) rows.
//!
//! E-values follow DIAMOND's ALP statistics ([`alp`]) rather than the plain
//! Karlin-Altschul formula, so they match the C++ program numerically.

pub mod alp;
pub mod alphabet;
pub mod cli;
pub mod db;
pub mod error;
pub mod gapped;
pub mod hsp;
pub mod io;
pub mod matrix;
pub mod matrix_data;
pub mod output;
pub mod reduction;
pub mod search;
pub mod seed;
pub mod shape;
pub mod translate;
pub mod ungapped;
pub mod workflow;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// `diamond-rust v0.1.0` banner used by `--version` and the stderr log.
pub fn version_string() -> String {
    format!("diamond-rust v{VERSION}")
}
