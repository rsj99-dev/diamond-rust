//! The `.dmnd`-style binary database.
//!
//! Layout (all integers little-endian, no padding):
//!
//! ```text
//! header   80 bytes  magic, version, flags, counts and section offsets
//! seq_offs (n+1) u64  offsets of every sequence in the letters section
//! def_offs (n+1) u64  offsets of every FASTA header in the defs section
//! letters  u8 * L     encoded amino acids, back to back
//! defs     bytes      newline-free FASTA headers, back to back
//! ```
//!
//! Like DIAMOND 2.x, the seed index is *not* stored: it is built in memory at
//! search time for the shapes the selected sensitivity needs, which keeps the
//! file small and lets `--shapes`/`--sensitive` change without rebuilding.

use crate::alphabet;
use crate::bail;
use crate::error::{Error, Result};
use crate::io;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

pub const MAGIC: &[u8; 8] = b"DRDMND01";
pub const VERSION: u32 = 1;
/// magic(8) + version(4) + flags(4) + 7 * u64.
const HEADER_LEN: usize = 72;
const FLAG_PROTEIN: u32 = 1;

#[derive(Debug, Clone, Copy)]
pub struct DbBuildStats {
    pub sequences: usize,
    pub letters: u64,
    pub seconds: f64,
}

/// Append `.dmnd` when the user did not type an extension (DIAMOND's
/// `--no-auto-append` default behaviour).
pub fn with_extension(path: &Path) -> PathBuf {
    if path.extension().is_some() {
        path.to_path_buf()
    } else {
        let mut s = path.as_os_str().to_os_string();
        s.push(".dmnd");
        PathBuf::from(s)
    }
}

/// Resolve a database path: as given, or with `.dmnd` appended.
pub fn resolve(path: &Path) -> Result<PathBuf> {
    if path.exists() {
        return Ok(path.to_path_buf());
    }
    let alt = with_extension(path);
    if alt.exists() {
        return Ok(alt);
    }
    bail!("Database file not found: {}", path.display())
}

pub fn build_database(input: &Path, output: &Path, threads: usize) -> Result<DbBuildStats> {
    let start = Instant::now();
    let records = io::read_fasta(input)?;
    if threads > 1 {
        // Encoding is trivially parallel; keep the sequential path for small
        // inputs and single-thread runs to avoid pool setup overhead.
        if records.len() > 4096 {
            use rayon::prelude::*;
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| Error::msg(format!("cannot start thread pool: {e}")))?;
            let encoded: Vec<Vec<u8>> = pool.install(|| {
                records
                    .par_iter()
                    .map(|r| alphabet::encode(&r.raw))
                    .collect()
            });
            return write_db(output, &records, &encoded, start);
        }
    }
    let encoded: Vec<Vec<u8>> = records.iter().map(|r| alphabet::encode(&r.raw)).collect();
    write_db(output, &records, &encoded, start)
}

fn write_db(
    output: &Path,
    records: &[io::Record],
    encoded: &[Vec<u8>],
    start: Instant,
) -> Result<DbBuildStats> {
    let n = records.len();
    if n == 0 {
        bail!("{}: no sequences", output.display());
    }
    let mut seq_offsets = Vec::with_capacity(n + 1);
    let mut def_offsets = Vec::with_capacity(n + 1);
    let mut letters: u64 = 0;
    let mut defs: u64 = 0;
    for (rec, seq) in records.iter().zip(encoded) {
        if seq.len() > u32::MAX as usize {
            bail!("sequence {} is too long for the database format", rec.id);
        }
        seq_offsets.push(letters);
        letters += seq.len() as u64;
        def_offsets.push(defs);
        defs += rec.header.len() as u64;
    }
    seq_offsets.push(letters);
    def_offsets.push(defs);

    let seq_off_start = HEADER_LEN as u64;
    let def_off_start = seq_off_start + 8 * (n as u64 + 1);
    let letters_start = def_off_start + 8 * (n as u64 + 1);
    let defs_start = letters_start + letters;

    let file = std::fs::File::create(output)
        .map_err(|e| Error::msg(format!("cannot create {}: {e}", output.display())))?;
    let mut out = std::io::BufWriter::with_capacity(1 << 20, file);

    let mut header = Vec::with_capacity(HEADER_LEN);
    header.extend_from_slice(MAGIC);
    header.extend_from_slice(&VERSION.to_le_bytes());
    header.extend_from_slice(&FLAG_PROTEIN.to_le_bytes());
    header.extend_from_slice(&(n as u64).to_le_bytes());
    header.extend_from_slice(&letters.to_le_bytes());
    header.extend_from_slice(&defs.to_le_bytes());
    header.extend_from_slice(&seq_off_start.to_le_bytes());
    header.extend_from_slice(&def_off_start.to_le_bytes());
    header.extend_from_slice(&letters_start.to_le_bytes());
    header.extend_from_slice(&defs_start.to_le_bytes());
    debug_assert_eq!(header.len(), HEADER_LEN);
    out.write_all(&header)?;

    for &off in &seq_offsets {
        out.write_all(&off.to_le_bytes())?;
    }
    for &off in &def_offsets {
        out.write_all(&off.to_le_bytes())?;
    }
    for seq in encoded {
        out.write_all(seq)?;
    }
    for rec in records {
        out.write_all(rec.header.as_bytes())?;
    }
    out.flush()?;
    drop(out);

    Ok(DbBuildStats {
        sequences: n,
        letters,
        seconds: start.elapsed().as_secs_f64(),
    })
}

/// A read-only, zero-copy view of a database file.
pub struct Database {
    data: Vec<u8>,
    path: PathBuf,
    seq_count: usize,
    total_letters: u64,
    seq_off_start: usize,
    def_off_start: usize,
    letters_start: usize,
    defs_start: usize,
    protein: bool,
}

impl Database {
    pub fn open(path: &Path) -> Result<Database> {
        let path = resolve(path)?;
        let data = std::fs::read(&path)
            .map_err(|e| Error::msg(format!("cannot read {}: {e}", path.display())))?;
        if data.len() < HEADER_LEN {
            bail!(
                "{}: not a diamond-rust database (file too short)",
                path.display()
            );
        }
        if &data[..8] != MAGIC {
            bail!(
                "{}: not a diamond-rust database (bad magic). Build one with `diamond-rust makedb`.",
                path.display()
            );
        }
        let version = u32::from_le_bytes(data[8..12].try_into().unwrap());
        if version != VERSION {
            bail!(
                "{}: database version {} is not supported (expected {})",
                path.display(),
                version,
                VERSION
            );
        }
        let flags = u32::from_le_bytes(data[12..16].try_into().unwrap());
        let seq_count = read_u64(&data, 16) as usize;
        let total_letters = read_u64(&data, 24);
        let defs_len = read_u64(&data, 32);
        let seq_off_start = read_u64(&data, 40) as usize;
        let def_off_start = read_u64(&data, 48) as usize;
        let letters_start = read_u64(&data, 56) as usize;
        let defs_start = read_u64(&data, 64) as usize;
        if defs_start as u64 + defs_len != data.len() as u64 {
            bail!("{}: truncated or corrupt database", path.display());
        }
        if seq_off_start + 8 * (seq_count + 1) > data.len()
            || def_off_start + 8 * (seq_count + 1) > data.len()
        {
            bail!("{}: corrupt offset table", path.display());
        }
        Ok(Database {
            data,
            path,
            seq_count,
            total_letters,
            seq_off_start,
            def_off_start,
            letters_start,
            defs_start,
            protein: flags & FLAG_PROTEIN != 0,
        })
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.seq_count
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.seq_count == 0
    }

    #[inline]
    pub fn total_letters(&self) -> u64 {
        self.total_letters
    }

    pub fn is_protein(&self) -> bool {
        self.protein
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    #[inline]
    fn seq_offset(&self, i: usize) -> usize {
        read_u64(&self.data, self.seq_off_start + 8 * i) as usize
    }

    #[inline]
    fn def_offset(&self, i: usize) -> usize {
        read_u64(&self.data, self.def_off_start + 8 * i) as usize
    }

    #[inline]
    pub fn seq(&self, i: usize) -> &[u8] {
        let lo = self.letters_start + self.seq_offset(i);
        let hi = self.letters_start + self.seq_offset(i + 1);
        &self.data[lo..hi]
    }

    #[inline]
    pub fn seq_len(&self, i: usize) -> usize {
        self.seq_offset(i + 1) - self.seq_offset(i)
    }

    #[inline]
    pub fn header(&self, i: usize) -> &str {
        let lo = self.defs_start + self.def_offset(i);
        let hi = self.defs_start + self.def_offset(i + 1);
        std::str::from_utf8(&self.data[lo..hi]).unwrap_or("")
    }

    #[inline]
    pub fn id(&self, i: usize) -> &str {
        alphabet::seq_id(self.header(i))
    }

    pub fn report(&self) -> String {
        format!(
            "Database: {}\nSequences: {}\nLetters: {}\nType: protein",
            self.path.display(),
            self.seq_count,
            fmt_count(self.total_letters)
        )
    }
}

#[inline]
fn read_u64(data: &[u8], off: usize) -> u64 {
    u64::from_le_bytes(data[off..off + 8].try_into().unwrap())
}

pub fn fmt_count(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_path(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!(
            "diamond_rust_test_{}_{}.dmnd",
            std::process::id(),
            name
        ));
        p
    }

    #[test]
    fn roundtrip_database() {
        let fasta = temp_path("in.fa");
        let db = temp_path("db");
        let mut f = std::fs::File::create(&fasta).unwrap();
        writeln!(f, ">s1 first").unwrap();
        writeln!(f, "MKWVTFISLL").unwrap();
        writeln!(f, ">s2 second").unwrap();
        writeln!(f, "ACDEFGHIKL").unwrap();
        drop(f);
        let stats = build_database(&fasta, &db, 1).unwrap();
        assert_eq!(stats.sequences, 2);
        assert_eq!(stats.letters, 20);
        let db = Database::open(&db).unwrap();
        assert_eq!(db.len(), 2);
        assert_eq!(db.id(0), "s1");
        assert_eq!(db.header(1), "s2 second");
        assert_eq!(db.seq(0), &alphabet::encode(b"MKWVTFISLL")[..]);
        assert_eq!(db.seq_len(1), 10);
        assert_eq!(db.total_letters(), 20);
        let _ = std::fs::remove_file(&fasta);
        let _ = std::fs::remove_file(&db.path);
    }
}
