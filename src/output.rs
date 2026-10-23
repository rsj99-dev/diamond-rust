//! BLAST tabular output (`--outfmt 6 [fields...]`).

use crate::bail;
use crate::db::Database;
use crate::error::{Error, Result};
use crate::hsp::Hsp;
use std::io::{BufWriter, Write};

/// The fields understood by `--outfmt 6 <field> ...`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Field {
    QSeqId,
    QLen,
    SSeqId,
    SLen,
    QStart,
    QEnd,
    SStart,
    SEnd,
    QSeq,
    SSeq,
    EValue,
    BitScore,
    Score,
    Length,
    PIdent,
    NIdent,
    Mismatch,
    Positive,
    GapOpen,
    Gaps,
    PPos,
    QFrame,
    QCovHsp,
    SCovHsp,
    Cigar,
    BTop,
    QTitle,
    STitle,
    FullSSeq,
    QSeqTranslated,
    QStrand,
}

impl Field {
    pub fn parse(name: &str) -> Option<Field> {
        Some(match name {
            "qseqid" => Field::QSeqId,
            "qlen" => Field::QLen,
            "sseqid" => Field::SSeqId,
            "slen" => Field::SLen,
            "qstart" => Field::QStart,
            "qend" => Field::QEnd,
            "sstart" => Field::SStart,
            "send" => Field::SEnd,
            "qseq" => Field::QSeq,
            "sseq" => Field::SSeq,
            "evalue" => Field::EValue,
            "bitscore" => Field::BitScore,
            "score" => Field::Score,
            "length" => Field::Length,
            "pident" => Field::PIdent,
            "nident" => Field::NIdent,
            "mismatch" => Field::Mismatch,
            "positive" => Field::Positive,
            "gapopen" => Field::GapOpen,
            "gaps" => Field::Gaps,
            "ppos" => Field::PPos,
            "qframe" => Field::QFrame,
            "qcovhsp" => Field::QCovHsp,
            "scovhsp" => Field::SCovHsp,
            "cigar" => Field::Cigar,
            "btop" => Field::BTop,
            "qtitle" => Field::QTitle,
            "stitle" => Field::STitle,
            "full_sseq" => Field::FullSSeq,
            "qseq_translated" => Field::QSeqTranslated,
            "qstrand" => Field::QStrand,
            _ => return None,
        })
    }

    /// The twelve columns of a bare `--outfmt 6`.
    pub fn defaults() -> Vec<Field> {
        vec![
            Field::QSeqId,
            Field::SSeqId,
            Field::PIdent,
            Field::Length,
            Field::Mismatch,
            Field::GapOpen,
            Field::QStart,
            Field::QEnd,
            Field::SStart,
            Field::SEnd,
            Field::EValue,
            Field::BitScore,
        ]
    }
}

/// Query metadata needed by the fields that are not part of the HSP.
pub struct QueryInfo<'a> {
    pub id: &'a str,
    pub header: &'a str,
    /// Encoded letters of the frame the HSP belongs to.
    pub letters: &'a [u8],
    /// Nucleotide length for `blastx`, amino-acid length for `blastp`.
    pub dna_len: usize,
    /// Length of the translated frame the HSP came from.
    pub frame_len: usize,
    pub translated: bool,
}

/// Parse the arguments of `--outfmt` (`6`, or `6 qseqid sseqid ...`).
pub fn parse_outfmt(tokens: &[String]) -> Result<Vec<Field>> {
    let Some(first) = tokens.first() else {
        return Ok(Field::defaults());
    };
    match first.as_str() {
        "6" => {}
        "0" | "5" | "102" => bail!("--outfmt {first} is not supported, use 6 (tabular)"),
        other => bail!("invalid --outfmt value {other:?}"),
    }
    if tokens.len() <= 1 {
        return Ok(Field::defaults());
    }
    let mut fields = Vec::with_capacity(tokens.len() - 1);
    for name in &tokens[1..] {
        let field = Field::parse(name)
            .ok_or_else(|| Error::msg(format!("Invalid output field: {name}")))?;
        fields.push(field);
    }
    Ok(fields)
}

pub struct Tabular {
    out: BufWriter<Box<dyn Write>>,
    fields: Vec<Field>,
}

impl Tabular {
    pub fn new(out: Box<dyn Write>, fields: Vec<Field>) -> Tabular {
        Tabular {
            out: BufWriter::with_capacity(1 << 20, out),
            fields,
        }
    }

    pub fn write(&mut self, h: &Hsp, q: &QueryInfo, db: &Database) -> Result<()> {
        let subject = db.seq(h.subject);
        let subject_len = db.seq_len(h.subject);
        let (qb, qe) = if q.translated {
            h.nt_range(q.dna_len)
        } else {
            (h.q_begin, h.q_end)
        };
        let query_len = if q.translated { q.dna_len } else { q.frame_len };

        let mut row: Vec<String> = Vec::with_capacity(self.fields.len());
        for &field in &self.fields {
            let value = match field {
                Field::QSeqId => q.id.to_string(),
                Field::QLen => query_len.to_string(),
                Field::SSeqId => db.id(h.subject).to_string(),
                Field::SLen => subject_len.to_string(),
                Field::QStart => (qb + 1).to_string(),
                Field::QEnd => qe.to_string(),
                Field::SStart => (h.s_begin + 1).to_string(),
                Field::SEnd => h.s_end.to_string(),
                Field::QSeq => h.aligned_strings(q.letters, subject).0,
                Field::SSeq => h.aligned_strings(q.letters, subject).1,
                Field::EValue => fmt_evalue(h.evalue),
                Field::BitScore => fmt_double(h.bits),
                Field::Score => h.score.to_string(),
                Field::Length => h.stats.length.to_string(),
                Field::PIdent => fmt_double(h.pident()),
                Field::NIdent => h.stats.identities.to_string(),
                Field::Mismatch => h.stats.mismatches.to_string(),
                Field::Positive => h.stats.positives.to_string(),
                Field::GapOpen => h.stats.gap_openings.to_string(),
                Field::Gaps => h.stats.gaps.to_string(),
                Field::PPos => fmt_double(h.ppos()),
                Field::QFrame => h.q_frame().to_string(),
                Field::QCovHsp => {
                    if query_len == 0 {
                        "0.0".to_string()
                    } else {
                        fmt_double((qe - qb) as f64 * 100.0 / query_len as f64)
                    }
                }
                Field::SCovHsp => {
                    if subject_len == 0 {
                        "0.0".to_string()
                    } else {
                        fmt_double((h.s_end - h.s_begin) as f64 * 100.0 / subject_len as f64)
                    }
                }
                Field::Cigar => h.cigar(),
                Field::BTop => h.btop(q.letters, subject),
                Field::QTitle => q.header.to_string(),
                Field::STitle => db.header(h.subject).to_string(),
                Field::FullSSeq => crate::alphabet::decode(subject),
                Field::QSeqTranslated => {
                    let letters = q.letters;
                    let end = h.q_end.min(letters.len());
                    let begin = h.q_begin.min(end);
                    crate::alphabet::decode(&letters[begin..end])
                }
                Field::QStrand => {
                    if h.is_reverse() {
                        "-".to_string()
                    } else {
                        "+".to_string()
                    }
                }
            };
            row.push(value);
        }
        self.out.write_all(row.join("\t").as_bytes())?;
        self.out.write_all(b"\n")?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.out.flush()?;
        Ok(())
    }
}

/// DIAMOND's `Util::String::format_double`.
///
/// Values of 100 or more are truncated to an integer (BLAST compatibility),
/// smaller ones keep one decimal.
pub fn fmt_double(x: f64) -> String {
    if x >= 100.0 {
        format!("{}", x.floor() as i64)
    } else {
        format!("{x:.1}")
    }
}

/// DIAMOND-style E-value formatting (`TextBuffer::print_e`, i.e. `%.2e`).
pub fn fmt_evalue(e: f64) -> String {
    if !e.is_finite() || e <= 0.0 {
        "0.0".to_string()
    } else {
        format!("{e:.2e}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_fields_match_blast() {
        let f = parse_outfmt(&["6".to_string()]).unwrap();
        assert_eq!(f.len(), 12);
        assert_eq!(f[0], Field::QSeqId);
        assert_eq!(f[11], Field::BitScore);
    }

    #[test]
    fn custom_fields_are_parsed() {
        let f = parse_outfmt(&[
            "6".to_string(),
            "qseqid".to_string(),
            "qframe".to_string(),
            "btop".to_string(),
        ])
        .unwrap();
        assert_eq!(f, vec![Field::QSeqId, Field::QFrame, Field::BTop]);
    }

    #[test]
    fn unknown_field_is_rejected() {
        assert!(parse_outfmt(&["6".to_string(), "nope".to_string()]).is_err());
        assert!(parse_outfmt(&["7".to_string()]).is_err());
    }

    #[test]
    fn evalues_are_formatted_like_diamond() {
        assert_eq!(fmt_evalue(0.0), "0.0");
        assert_eq!(fmt_evalue(0.001), "1.00e-3");
        assert_eq!(fmt_evalue(4.7e-45), "4.70e-45");
        assert_eq!(fmt_evalue(3.21e-105), "3.21e-105");
    }

    #[test]
    fn doubles_are_formatted_like_diamond() {
        assert_eq!(fmt_double(319.7), "319");
        assert_eq!(fmt_double(105.49), "105");
        assert_eq!(fmt_double(80.87), "80.9");
        assert_eq!(fmt_double(22.581), "22.6");
    }
}
