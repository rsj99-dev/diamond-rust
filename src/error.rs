//! Error type and `bail!` helper.
//!
//! Every fallible path in the crate reports a message that is already
//! user-facing, so the CLI can print it verbatim (DIAMOND style).

use std::fmt;
use std::io;
use std::num::{ParseFloatError, ParseIntError};

#[derive(Debug, Clone)]
pub struct Error(String);

impl Error {
    pub fn msg(m: impl Into<String>) -> Error {
        Error(m.into())
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Error {
        Error(e.to_string())
    }
}

impl From<ParseIntError> for Error {
    fn from(e: ParseIntError) -> Error {
        Error(e.to_string())
    }
}

impl From<ParseFloatError> for Error {
    fn from(e: ParseFloatError) -> Error {
        Error(e.to_string())
    }
}

pub type Result<T> = std::result::Result<T, Error>;

#[macro_export]
macro_rules! bail {
    ($($arg:tt)*) => {
        return Err($crate::error::Error::msg(format!($($arg)*)))
    };
}
