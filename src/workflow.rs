//! Command dispatch and the end-to-end workflows.

use crate::alphabet;
use crate::cli;
use crate::db::{self, Database};
use crate::error::{Error, Result};
use crate::io;
use crate::matrix::ScoreMatrix;
use crate::output::{self, QueryInfo, Tabular};
use crate::search::{self, SearchConfig, Target};
use crate::translate;
use std::io::Write;
use std::time::Instant;

pub fn run() -> Result<i32> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(command) = args.first() else {
        print_usage();
        return Ok(1);
    };
    match command.as_str() {
        "version" | "--version" => {
            println!("{}", crate::version_string());
            Ok(0)
        }
        "help" | "--help" | "-h" => {
            print_usage();
            Ok(0)
        }
        "makedb" => cmd_makedb(&args[1..]),
        "blastp" => cmd_blast(&args[1..], false),
        "blastx" => cmd_blast(&args[1..], true),
        "dbinfo" => cmd_dbinfo(&args[1..]),
        other => {
            eprintln!("Unknown command: {other}\n");
            print_usage();
            Ok(1)
        }
    }
}

fn cmd_makedb(args: &[String]) -> Result<i32> {
    let opts = cli::parse_makedb(args)?;
    let output = if opts.no_auto_append {
        opts.output.clone()
    } else {
        db::with_extension(&opts.output)
    };
    if !opts.quiet {
        eprintln!("{}", crate::version_string());
    }
    let stats = db::build_database(&opts.input, &output, opts.threads)?;
    if !opts.quiet {
        eprintln!("Database: {}", output.display());
        eprintln!("Sequences: {}", stats.sequences);
        eprintln!("Letters: {}", db::fmt_count(stats.letters));
        eprintln!("Total time = {:.3}s", stats.seconds);
    }
    Ok(0)
}

fn cmd_dbinfo(args: &[String]) -> Result<i32> {
    let mut db_path: Option<std::path::PathBuf> = None;
    let mut i = 0usize;
    while i < args.len() {
        let arg = &args[i];
        let (name, inline) = match arg.strip_prefix("--").or_else(|| arg.strip_prefix('-')) {
            Some(rest) => match rest.split_once('=') {
                Some((n, v)) => (n, Some(v)),
                None => (rest, None),
            },
            None => return Err(Error::msg(format!("Unexpected argument: {arg}"))),
        };
        match name {
            "db" | "d" => {
                let v = match inline {
                    Some(v) => v.to_string(),
                    None => {
                        i += 1;
                        args.get(i)
                            .cloned()
                            .ok_or_else(|| Error::msg("Option --db requires a value"))?
                    }
                };
                db_path = Some(std::path::PathBuf::from(v));
            }
            other => return Err(Error::msg(format!("Unknown option: --{other}"))),
        }
        i += 1;
    }
    let path = db_path.ok_or_else(|| Error::msg("Missing option: --db"))?;
    let db = Database::open(&path)?;
    println!("{}", db.report());
    Ok(0)
}

fn cmd_blast(args: &[String], translated: bool) -> Result<i32> {
    let opts = cli::parse_blast(args, translated)?;
    if opts.verbose && !opts.quiet {
        eprintln!("{}", crate::version_string());
    }
    let start = Instant::now();

    let database = Database::open(&opts.db)?;
    let records = io::read_fasta(&opts.query)?;
    let matrix = ScoreMatrix::standard(&opts.matrix, opts.gap_open, opts.gap_extend)?;

    // Translations for blastx, plain encoded residues for blastp.
    let mut targets: Vec<Target> = Vec::new();
    let mut target_of: Vec<Vec<Option<usize>>> = vec![vec![None; 6]; records.len()];
    for (qi, record) in records.iter().enumerate() {
        if translated {
            let min_orf = translate::min_orf_len(record.raw.len(), opts.min_orf);
            for (frame_index, mut letters) in
                translate::six_frames(&record.raw).into_iter().enumerate()
            {
                translate::mask_short_orfs(&mut letters, min_orf);
                if letters.iter().filter(|&&l| l < alphabet::X_LETTER).count() == 0 {
                    continue;
                }
                target_of[qi][frame_index] = Some(targets.len());
                targets.push(Target {
                    query: qi,
                    frame: frame_index as u8,
                    letters,
                });
            }
        } else {
            let letters = alphabet::encode(&record.raw);
            if letters.is_empty() {
                continue;
            }
            target_of[qi][0] = Some(targets.len());
            targets.push(Target {
                query: qi,
                frame: 0,
                letters,
            });
        }
    }

    let mut cfg = SearchConfig::new(matrix);
    cfg.evalue = opts.evalue;
    cfg.max_target_seqs = opts.max_target_seqs;
    cfg.max_hsps = opts.max_hsps;
    cfg.sensitivity = opts.sensitivity.clone();
    cfg.threads = opts.threads;
    cfg.index_chunks = opts.index_chunks;
    cfg.seed_cut = opts.seed_cut;

    let results = search::run(&database, &targets, &cfg)?;
    let per_query = search::merge_frames(
        results,
        &targets,
        records.len(),
        cfg.max_hsps,
        cfg.max_target_seqs,
    );

    let fields = output::parse_outfmt(&opts.outfmt)?;
    let sink: Box<dyn Write> = match &opts.out {
        Some(path) => Box::new(
            std::fs::File::create(path)
                .map_err(|e| Error::msg(format!("cannot create {}: {e}", path.display())))?,
        ),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = Tabular::new(sink, fields);
    let mut reported = 0usize;
    let mut queries_hit = 0usize;
    for (qi, hsps) in per_query.iter().enumerate() {
        if hsps.is_empty() {
            continue;
        }
        queries_hit += 1;
        let record = &records[qi];
        for hsp in hsps {
            let Some(Some(ti)) = target_of.get(qi).map(|v| v[hsp.frame as usize]) else {
                continue;
            };
            let letters = &targets[ti].letters;
            let info = QueryInfo {
                id: &record.id,
                header: &record.header,
                letters,
                dna_len: record.raw.len(),
                frame_len: letters.len(),
                translated,
            };
            writer.write(hsp, &info, &database)?;
            reported += 1;
        }
    }
    writer.flush()?;

    if !opts.quiet {
        eprintln!("Total time = {:.3}s", start.elapsed().as_secs_f64());
        eprintln!(
            "Reported {} pairwise alignments, {} HSPs.",
            reported, reported
        );
        eprintln!("{} queries aligned.", queries_hit);
    }
    Ok(0)
}

fn print_usage() {
    println!("{}", usage_text());
}

fn usage_text() -> String {
    format!(
        "{} — a Rust reimplementation of DIAMOND\n\
         \n\
         Usage:\n\
         \x20 diamond-rust makedb --in <FASTA> --db <DB> [options]\n\
         \x20 diamond-rust blastp --db <DB> --query <FASTA> [options]\n\
         \x20 diamond-rust blastx --db <DB> --query <FASTA> [options]\n\
         \x20 diamond-rust dbinfo --db <DB>\n\
         \n\
         Options:\n\
         \x20 -o, --out <FILE>          output file (default: stdout)\n\
         \x20 -f, --outfmt <FORMAT>     output format: 6 [field ...] (default: 6)\n\
         \x20 -e, --evalue <VALUE>      maximum E-value to report (default: 0.001)\n\
         \x20 -k, --max-target-seqs <N> maximum number of target sequences (default: 25)\n\
         \x20     --max-hsps <N>        maximum number of HSPs per target (default: 1)\n\
         \x20 -p, --threads <N>         number of CPU threads (default: all cores)\n\
         \x20     --index-chunks <N>    database chunks indexed at a time (default: preset specific)\n\
         \x20     --matrix <NAME>       scoring matrix (default: BLOSUM62)\n\
         \x20     --gapopen <N>         gap open penalty (default: matrix specific)\n\
         \x20     --gapextend <N>       gap extension penalty (default: matrix specific)\n\
         \x20 -l, --min-orf <N>         minimum ORF length for blastx\n\
         \x20     --seed-cut <VALUE>    seed complexity cutoff (default: preset specific)\n\
         \x20     --fast | --faster | --mid-sensitive | --sensitive\n\
         \x20     --more-sensitive | --very-sensitive | --ultra-sensitive\n\
         \x20 -q, --quiet              suppress the progress log\n\
         \n\
         blastx additionally accepts --query-gencode (translation tables other\n\
         than the standard code are not implemented).",
        crate::version_string()
    )
}

#[cfg(test)]
mod tests {
    use super::usage_text;

    #[test]
    fn usage_text_mentions_every_command() {
        // Keep the help text in sync with the dispatcher.
        let text = usage_text();
        for cmd in ["makedb", "blastp", "blastx", "dbinfo"] {
            assert!(text.contains(cmd), "usage text is missing `{cmd}`");
        }
    }
}
