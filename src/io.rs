//! FASTA reading.
//!
//! Records keep their raw (upper-cased, layout-stripped) bytes; alphabet
//! encoding happens where it is needed (database build, query loading,
//! six-frame translation).

use crate::error::{Error, Result};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct Record {
    /// Full FASTA header without the leading `>`.
    pub header: String,
    /// First whitespace-delimited token of the header.
    pub id: String,
    /// Upper-case sequence with whitespace, `-` and `.` removed.
    pub raw: Vec<u8>,
}

impl Record {
    pub fn len(&self) -> usize {
        self.raw.len()
    }

    pub fn is_empty(&self) -> bool {
        self.raw.is_empty()
    }
}

pub fn read_fasta(path: &Path) -> Result<Vec<Record>> {
    let data = std::fs::read(path)
        .map_err(|e| Error::msg(format!("cannot read {}: {e}", path.display())))?;
    let records = parse_fasta(&data)?;
    if records.is_empty() {
        crate::bail!("{}: no FASTA records found", path.display());
    }
    Ok(records)
}

pub fn parse_fasta(data: &[u8]) -> Result<Vec<Record>> {
    let mut out: Vec<Record> = Vec::new();
    let mut header: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    let mut seen_header = false;

    for line in split_lines(data) {
        if line.first() == Some(&b'>') {
            if let Some(h) = header.take() {
                out.push(finish(h, &mut seq));
            }
            seen_header = true;
            let text = String::from_utf8_lossy(&line[1..]).trim().to_string();
            header = Some(text);
        } else if header.is_some() {
            for &b in line {
                if b.is_ascii_whitespace() || b == b'-' || b == b'.' {
                    continue;
                }
                seq.push(b.to_ascii_uppercase());
            }
        } else if !line.is_empty() && line[0] != b';' && line[0] != b'#' {
            return Err(Error::msg(
                "malformed FASTA: sequence data found before the first '>' header",
            ));
        }
    }
    if let Some(h) = header.take() {
        out.push(finish(h, &mut seq));
    }
    if !seen_header && !data.is_empty() {
        crate::bail!("malformed FASTA: no '>' header found");
    }
    Ok(out)
}

fn finish(header: String, seq: &mut Vec<u8>) -> Record {
    let id = header
        .split(|c: char| c.is_ascii_whitespace())
        .next()
        .unwrap_or("")
        .to_string();
    let id = if id.is_empty() {
        format!("seq{}", seq.len().max(1))
    } else {
        id
    };
    Record {
        header,
        id,
        raw: std::mem::take(seq),
    }
}

fn split_lines(data: &[u8]) -> impl Iterator<Item = &[u8]> {
    data.split(|&b| b == b'\n').map(|line| {
        if line.last() == Some(&b'\r') {
            &line[..line.len() - 1]
        } else {
            line
        }
    })
}

/// Heuristic alphabet guess, mirrored from `makeblastdb -dbtype auto`.
/// A sequence is treated as nucleotide when >= 85% of its sampled letters are
/// A/C/G/T/U/N.
pub fn looks_like_dna(records: &[Record]) -> bool {
    let mut canonical = 0usize;
    let mut total = 0usize;
    for rec in records {
        for &b in rec.raw.iter().take(20_000) {
            total += 1;
            if matches!(b, b'A' | b'C' | b'G' | b'T' | b'U' | b'N') {
                canonical += 1;
            }
        }
        if total >= 20_000 {
            break;
        }
    }
    total > 0 && canonical * 100 >= total * 85
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multiline_records() {
        let recs =
            parse_fasta(b">seq1 first description\r\nACGT\r\nACGT\r\n>seq2\nAAAA\n").unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "seq1");
        assert_eq!(recs[0].header, "seq1 first description");
        assert_eq!(recs[0].raw, b"ACGTACGT");
        assert_eq!(recs[1].id, "seq2");
        assert_eq!(recs[1].raw, b"AAAA");
    }

    #[test]
    fn rejects_sequence_before_header() {
        assert!(parse_fasta(b"ACGT\n>a\nAC\n").is_err());
    }

    #[test]
    fn detects_dna() {
        let dna = parse_fasta(b">a\nACGTACGTTT\n").unwrap();
        assert!(looks_like_dna(&dna));
        let prot = parse_fasta(b">a\nMKWVTFISLLFLFSSAYSR\n").unwrap();
        assert!(!looks_like_dna(&prot));
    }
}
