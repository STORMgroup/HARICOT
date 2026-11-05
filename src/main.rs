//! Binary entrypoint. All logic lives in the library (`haricot::cli`) so it
//! stays unit-testable.

fn main() -> anyhow::Result<()> {
    haricot::cli::run()
}
