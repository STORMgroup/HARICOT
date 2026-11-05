//! The strict command line of the pipeline programs that read their flags by name (`haricot-pair-features`,
//! `haricot-graph-features`, `haricot-forest`, `haricot-rescue`; `#[path]`-included like `cache.rs`). A program declares its
//! flags, those that take a value and the switches; an undeclared `--` flag, a flag without its value,
//! a flag given twice and a value that does not parse are errors, as in `haricot-detect` and `haricot-layout`:
//! a typo stops the program instead of leaving a default in place. Every other argument is
//! positional, in order.
#![allow(dead_code)]

use std::collections::HashMap;
use std::fmt::Display;
use std::str::FromStr;

use anyhow::{Context, Result, anyhow, bail, ensure};

/// A program's parsed command line.
pub struct Cli {
    /// The positional arguments, in order.
    pub pos: Vec<String>,
    vals: HashMap<&'static str, String>,
    on: Vec<&'static str>,
    valued: &'static [&'static str],
    switches: &'static [&'static str],
}

impl Cli {
    /// Parse `argv` (the arguments after the program's name, or after its subcommand): `valued` are
    /// the flags that take the next argument as their value, `switches` the flags that take none.
    pub fn parse(
        argv: &[String],
        valued: &'static [&'static str],
        switches: &'static [&'static str],
        usage: &str,
    ) -> Result<Cli> {
        let mut c = Cli {
            pos: Vec::new(),
            vals: HashMap::new(),
            on: Vec::new(),
            valued,
            switches,
        };
        let mut it = argv.iter();
        while let Some(a) = it.next() {
            if let Some(&n) = valued.iter().find(|&&n| n == a) {
                let v = it
                    .next()
                    .with_context(|| format!("{n} needs a value; {usage}"))?;
                ensure!(
                    c.vals.insert(n, v.clone()).is_none(),
                    "{n} is given twice; {usage}"
                );
            } else if let Some(&n) = switches.iter().find(|&&n| n == a) {
                ensure!(!c.on.contains(&n), "{n} is given twice; {usage}");
                c.on.push(n);
            } else if a.starts_with("--") {
                bail!("unknown flag {a}; {usage}");
            } else {
                c.pos.push(a.clone());
            }
        }
        Ok(c)
    }

    /// Exactly `n` positional arguments.
    pub fn positional(&self, n: usize, usage: &str) -> Result<()> {
        ensure!(
            self.pos.len() == n,
            "{} positional arguments ({:?}), not {n}; {usage}",
            self.pos.len(),
            self.pos
        );
        Ok(())
    }

    /// The value of the flag `n`, if given (`n` must be declared as taking one).
    pub fn opt(&self, n: &str) -> Option<String> {
        assert!(self.valued.contains(&n), "{n}: not a declared flag");
        self.vals.get(n).cloned()
    }

    /// The value of the flag `n` parsed as a `T`, `d` when it is not given; a value that does not
    /// parse is an error.
    pub fn num<T: FromStr>(&self, n: &str, d: T) -> Result<T>
    where
        T::Err: Display,
    {
        match self.opt(n) {
            None => Ok(d),
            Some(v) => v.parse().map_err(|e| anyhow!("{n} {v}: {e}")),
        }
    }

    /// Whether the switch `n` is given (`n` must be declared as a switch).
    pub fn has(&self, n: &str) -> bool {
        assert!(self.switches.contains(&n), "{n}: not a declared switch");
        self.on.contains(&n)
    }

    /// The canonical pair frame (`haricot::readorder::query_first`), the default:
    /// off with `--no-symmetric` (the owner's frame); `--symmetric` states the
    /// default. Both together are an error.
    pub fn symmetric(&self) -> Result<bool> {
        let (on, off) = (self.has("--symmetric"), self.has("--no-symmetric"));
        ensure!(!(on && off), "--symmetric and --no-symmetric together");
        Ok(!off)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const V: &[&str] = &["--out", "--band"];
    const S: &[&str] = &["--symmetric", "--no-symmetric", "--no-dtw"];

    fn p(s: &str) -> Result<Cli> {
        let v: Vec<String> = s.split_whitespace().map(String::from).collect();
        Cli::parse(&v, V, S, "usage")
    }

    #[test]
    fn flags_values_and_positionals() {
        let c = p("a.pod5 c.bin --out f.tsv --band 32 --no-dtw").unwrap();
        assert_eq!(c.pos, ["a.pod5", "c.bin"]);
        assert!(c.positional(2, "u").is_ok() && c.positional(3, "u").is_err());
        assert_eq!(c.opt("--out").as_deref(), Some("f.tsv"));
        assert_eq!(c.num("--band", 64i64).unwrap(), 32);
        assert!(c.has("--no-dtw") && !c.has("--symmetric"));
        assert!(c.symmetric().unwrap());
        assert_eq!(p("x").unwrap().num("--band", 64i64).unwrap(), 64);
    }

    #[test]
    fn typos_and_malformed_values_are_errors() {
        assert!(p("x --bnad 32").is_err());
        assert!(p("x --out").is_err());
        assert!(p("x --out a --out b").is_err());
        assert!(p("x --no-dtw --no-dtw").is_err());
        assert!(p("x --band 3x").unwrap().num("--band", 64i64).is_err());
        assert!(p("x --band 3.5").unwrap().num("--band", 64i64).is_err());
    }

    #[test]
    fn the_frame_is_canonical_unless_no_symmetric() {
        assert!(p("x --symmetric").unwrap().symmetric().unwrap());
        assert!(!p("x --no-symmetric").unwrap().symmetric().unwrap());
        assert!(
            p("x --symmetric --no-symmetric")
                .unwrap()
                .symmetric()
                .is_err()
        );
    }
}
