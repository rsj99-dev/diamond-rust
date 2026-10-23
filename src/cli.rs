//! Command line parsing.
//!
//! DIAMOND's option names (and their short forms) are accepted so existing
//! scripts keep working. Options that only tune internals we do not implement
//! (block size, temporary directory, composition-based statistics, ...) are
//! parsed and ignored.

use crate::bail;
use crate::error::{Error, Result};
use std::path::PathBuf;

pub fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

pub struct MakedbArgs {
    pub input: PathBuf,
    pub output: PathBuf,
    pub threads: usize,
    pub no_auto_append: bool,
    pub quiet: bool,
}

pub struct BlastArgs {
    pub db: PathBuf,
    pub query: PathBuf,
    pub out: Option<PathBuf>,
    pub outfmt: Vec<String>,
    pub evalue: f64,
    pub max_target_seqs: usize,
    pub max_hsps: usize,
    pub threads: usize,
    /// `--index-chunks`; `None` keeps the sensitivity preset's value.
    pub index_chunks: Option<usize>,
    pub sensitivity: String,
    pub matrix: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub min_orf: Option<usize>,
    pub seed_cut: Option<f64>,
    pub quiet: bool,
    pub verbose: bool,
    pub translated: bool,
}

/// Split `--name=value`, `--name`, `-n` into a name and an optional value.
fn split_opt(arg: &str) -> Option<(&str, Option<&str>)> {
    if !arg.starts_with('-') || arg == "-" {
        return None;
    }
    let s = arg.trim_start_matches('-');
    if s.is_empty() {
        return None;
    }
    match s.split_once('=') {
        Some((name, value)) => Some((name, Some(value))),
        None => Some((s, None)),
    }
}

fn value(args: &[String], i: &mut usize, inline: Option<&str>, name: &str) -> Result<String> {
    if let Some(v) = inline {
        return Ok(v.to_string());
    }
    *i += 1;
    args.get(*i)
        .cloned()
        .ok_or_else(|| Error::msg(format!("Option --{name} requires a value")))
}

fn parse_usize(text: &str, name: &str) -> Result<usize> {
    text.parse::<usize>()
        .map_err(|_| Error::msg(format!("Invalid value for --{name}: {text:?}")))
}

fn parse_i32(text: &str, name: &str) -> Result<i32> {
    text.parse::<i32>()
        .map_err(|_| Error::msg(format!("Invalid value for --{name}: {text:?}")))
}

fn parse_f64(text: &str, name: &str) -> Result<f64> {
    text.parse::<f64>()
        .map_err(|_| Error::msg(format!("Invalid value for --{name}: {text:?}")))
}

pub fn parse_makedb(args: &[String]) -> Result<MakedbArgs> {
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut threads = default_threads();
    let mut no_auto_append = false;
    let mut quiet = false;
    let mut i = 0usize;
    while i < args.len() {
        let arg = &args[i];
        let Some((name, inline)) = split_opt(arg) else {
            bail!("Unexpected argument: {arg}");
        };
        match name {
            "in" | "i" => input = Some(PathBuf::from(value(args, &mut i, inline, name)?)),
            "db" | "d" => output = Some(PathBuf::from(value(args, &mut i, inline, name)?)),
            "threads" | "p" => threads = parse_usize(&value(args, &mut i, inline, name)?, name)?,
            "no-auto-append" => no_auto_append = true,
            "quiet" => quiet = true,
            "verbose" => {}
            "compress" | "compression" => {
                let _ = value(args, &mut i, inline, name)?;
            }
            other => bail!("Unknown option: --{other}"),
        }
        i += 1;
    }
    let input = input.ok_or_else(|| Error::msg("Missing option: --in"))?;
    let output = output.ok_or_else(|| Error::msg("Missing option: --db"))?;
    Ok(MakedbArgs {
        input,
        output,
        threads: threads.max(1),
        no_auto_append,
        quiet,
    })
}

pub fn parse_blast(args: &[String], translated: bool) -> Result<BlastArgs> {
    let mut db: Option<PathBuf> = None;
    let mut query: Option<PathBuf> = None;
    let mut out: Option<PathBuf> = None;
    let mut outfmt: Vec<String> = Vec::new();
    let mut evalue = 0.001f64;
    let mut max_target_seqs = 25usize;
    let mut max_hsps = 1usize;
    let mut threads = default_threads();
    let mut index_chunks: Option<usize> = None;
    let mut sensitivity = "default".to_string();
    let mut matrix = "BLOSUM62".to_string();
    let mut gap_open = -1i32;
    let mut gap_extend = -1i32;
    let mut min_orf: Option<usize> = None;
    let mut seed_cut: Option<f64> = None;
    let mut quiet = false;
    let mut verbose = true;

    let mut i = 0usize;
    while i < args.len() {
        let arg = &args[i];
        let Some((name, inline)) = split_opt(arg) else {
            bail!("Unexpected argument: {arg}");
        };
        match name {
            "db" | "d" => db = Some(PathBuf::from(value(args, &mut i, inline, name)?)),
            "query" | "q" => query = Some(PathBuf::from(value(args, &mut i, inline, name)?)),
            "out" | "o" => out = Some(PathBuf::from(value(args, &mut i, inline, name)?)),
            "outfmt" | "f" => {
                outfmt.push(value(args, &mut i, inline, name)?);
                // `--outfmt 6 qseqid sseqid ...`: the field names follow the
                // format number as separate tokens.
                while i + 1 < args.len() && !args[i + 1].starts_with('-') {
                    i += 1;
                    outfmt.push(args[i].clone());
                }
            }
            "evalue" | "e" => evalue = parse_f64(&value(args, &mut i, inline, name)?, name)?,
            "max-target-seqs" | "k" => {
                max_target_seqs = parse_usize(&value(args, &mut i, inline, name)?, name)?
            }
            "max-hsps" => max_hsps = parse_usize(&value(args, &mut i, inline, name)?, name)?,
            "threads" | "p" => threads = parse_usize(&value(args, &mut i, inline, name)?, name)?,
            "index-chunks" => {
                index_chunks = Some(parse_usize(&value(args, &mut i, inline, name)?, name)?)
            }
            "matrix" => matrix = value(args, &mut i, inline, name)?,
            "gapopen" => gap_open = parse_i32(&value(args, &mut i, inline, name)?, name)?,
            "gapextend" => gap_extend = parse_i32(&value(args, &mut i, inline, name)?, name)?,
            "min-orf" | "l" => {
                min_orf = Some(parse_usize(&value(args, &mut i, inline, name)?, name)?)
            }
            "seed-cut" => seed_cut = Some(parse_f64(&value(args, &mut i, inline, name)?, name)?),
            "fast" => sensitivity = "fast".to_string(),
            "faster" => sensitivity = "faster".to_string(),
            "mid-sensitive" => sensitivity = "mid-sensitive".to_string(),
            "sensitive" => sensitivity = "sensitive".to_string(),
            "more-sensitive" => sensitivity = "more-sensitive".to_string(),
            "very-sensitive" => sensitivity = "very-sensitive".to_string(),
            "ultra-sensitive" => sensitivity = "ultra-sensitive".to_string(),
            "quiet" => quiet = true,
            "verbose" => verbose = true,
            // Options that only affect internals this implementation does not
            // model; accept and ignore them so DIAMOND command lines run.
            "query-gencode"
            | "block-size"
            | "tmpdir"
            | "comp-based-stats"
            | "masking"
            | "motif-masking"
            | "soft-masking"
            | "algo"
            | "bin"
            | "swipe"
            | "ext"
            | "query-cover"
            | "subject-cover"
            | "id"
            | "query-length"
            | "approx-id"
            | "log"
            | "compress"
            | "chunk-size"
            | "load-threads"
            | "minichunk"
            | "freq-sd"
            | "sketch-size"
            | "tile-size"
            | "id2"
            | "padding"
            | "narrow-band-cov"
            | "narrow-band-factor"
            | "min-query-len"
            | "seed-freq"
            | "matrix-scan"
            | "max-hsp"
            | "xdrop"
            | "gapped-xdrop"
            | "ungapped-evalue"
            | "gapped-filter-evalue"
            | "toppercent"
            | "top"
            | "sequence-split"
            | "blastx"
            | "frameshift"
            | "xml"
            | "json" => {
                if inline.is_none() && i + 1 < args.len() && !args[i + 1].starts_with('-') {
                    // Most of these take a value; consume it.
                    if !matches!(name, "swipe" | "xml" | "json") {
                        i += 1;
                    }
                }
            }
            "no-unmask" | "unal" | "header" | "no-reorder" | "ignore-warnings" | "self"
            | "iterate" | "no-auto-append" => {}
            other => bail!("Unknown option: --{other}"),
        }
        i += 1;
    }

    let db = db.ok_or_else(|| Error::msg("Missing option: --db"))?;
    let query = query.ok_or_else(|| Error::msg("Missing option: --query"))?;
    Ok(BlastArgs {
        db,
        query,
        out,
        outfmt,
        evalue,
        max_target_seqs: max_target_seqs.max(1),
        max_hsps: max_hsps.max(1),
        threads: threads.max(1),
        index_chunks,
        sensitivity,
        matrix,
        gap_open,
        gap_extend,
        min_orf,
        seed_cut,
        quiet,
        verbose,
        translated,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(s: &str) -> Vec<String> {
        s.split_whitespace().map(|x| x.to_string()).collect()
    }

    #[test]
    fn makedb_options() {
        let a = parse_makedb(&args("-i db.fa -d out -p 4 --no-auto-append")).unwrap();
        assert_eq!(a.input, PathBuf::from("db.fa"));
        assert_eq!(a.output, PathBuf::from("out"));
        assert_eq!(a.threads, 4);
        assert!(a.no_auto_append);
    }

    #[test]
    fn blast_options() {
        let a = parse_blast(
            &args(
                "--db nr --query q.fa -e 1e-5 -k 10 -p 8 --more-sensitive --outfmt 6 qseqid pident",
            ),
            false,
        )
        .unwrap();
        assert_eq!(a.evalue, 1e-5);
        assert_eq!(a.max_target_seqs, 10);
        assert_eq!(a.threads, 8);
        assert_eq!(a.sensitivity, "more-sensitive");
        assert_eq!(a.outfmt, vec!["6", "qseqid", "pident"]);
    }

    #[test]
    fn inline_values_are_supported() {
        let a = parse_blast(&args("--db=nr --query=q.fa --threads=2"), true).unwrap();
        assert_eq!(a.threads, 2);
        assert!(a.translated);
    }

    #[test]
    fn seed_cut_is_optional() {
        let a = parse_blast(&args("--db nr --query q.fa --seed-cut 0.55"), false).unwrap();
        assert_eq!(a.seed_cut, Some(0.55));
        let b = parse_blast(&args("--db nr --query q.fa"), false).unwrap();
        assert_eq!(b.seed_cut, None);
    }

    #[test]
    fn missing_required_options_are_errors() {
        assert!(parse_blast(&args("-q q.fa"), false).is_err());
        assert!(parse_blast(&args("-d nr"), false).is_err());
        assert!(parse_makedb(&args("-i x.fa")).is_err());
    }

    #[test]
    fn unknown_options_are_rejected() {
        assert!(parse_blast(&args("-d nr -q q.fa --nope"), false).is_err());
    }
}
