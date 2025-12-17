//! Spawn multiple threads, each with its own `Aligner`.
//!
//! `Aligner` and `Parameters` are intentionally `!Send`/`!Sync` in this crate, so you can't
//! create an aligner on one thread and move it to another. The usual pattern is to construct
//! one aligner inside each worker thread and keep it there for the thread's lifetime.
//!
//! You can verify the Rust wrapper does not trigger upstream abPOA global-table data races
//! under TSAN, e.g.:
//! ```
//! MallocNanoZone=0 TSAN_OPTIONS="halt_on_error=1" CC=clang CFLAGS_aarch64_apple_darwin="-fsanitize=thread -fno-omit-frame-pointer -O1 -g" RUSTFLAGS="-Zsanitizer=thread" cargo +nightly run -Zbuild-std --target aarch64-apple-darwin -p abpoa --example multi_thread_aligners
//! ```
//! This should run without ThreadSanitizer reporting a data race.

use abpoa::{Aligner, Parameters, Result, SequenceBatch};
use std::{sync::Arc, thread, time::Instant};

#[path = "shared/sim.rs"]
mod sim;

fn main() -> Result<()> {
    let threads = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let iters = 10;

    let seed = Some(42);

    let consensus = sim::make_consensus(sim::ConsensusOptions {
        len: 200,
        seed,
        ..Default::default()
    });
    let reads: Arc<Vec<Vec<u8>>> = Arc::new(
        sim::simulate_reads(
            &consensus,
            40,
            sim::ErrorRates::new(0.05, 0.0025, 0.0025),
            seed,
        )
        .collect::<Vec<Vec<u8>>>(),
    );

    let expected = run_once(iters, reads.clone(), None, false)?;
    let mut handles = Vec::with_capacity(threads);

    eprintln!("spawning {threads} worker threads x {iters} iterations");

    for worker_idx in 0..threads {
        let reads = reads.clone();
        handles.push(thread::spawn(move || {
            run_once(iters, reads, Some(worker_idx), true)
        }));
    }

    for (idx, handle) in handles.into_iter().enumerate() {
        let got = handle.join().expect("worker thread panicked")?;
        if got != expected {
            return Err(abpoa::Error::InvalidInput(
                format!("thread {idx} consensus mismatch: expected={expected:?} got={got:?}")
                    .into(),
            ));
        }
    }

    println!(
        "ok: {threads} threads x {iters} iterations, consensus_len={}",
        expected.len()
    );
    Ok(())
}

fn run_once(
    iters: usize,
    reads: Arc<Vec<Vec<u8>>>,
    worker_idx: Option<usize>,
    verbose: bool,
) -> Result<String> {
    let read_refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
    let mut params = Parameters::configure()?;
    params.set_consensus(abpoa::ConsensusAlgorithm::MostFrequent, 2, 0.3)?;
    let mut aligner = Aligner::with_params(params)?;

    let label = match worker_idx {
        Some(idx) => format!("worker-{idx}"),
        None => "main".to_string(),
    };
    let thread_id = thread::current().id();
    let start = Instant::now();

    if verbose {
        eprintln!(
            "[{label}] start: thread_id={thread_id:?}, reads={}",
            read_refs.len()
        );
    }

    let mut expected: Option<String> = None;
    for iter_idx in 0..iters {
        let iter_start = Instant::now();
        let result = aligner.msa(SequenceBatch::from_sequences(&read_refs))?;
        let got = result
            .clusters
            .first()
            .map(|c| c.consensus.clone())
            .unwrap_or_default();

        if verbose {
            eprintln!(
                "[{label}] iter {}/{}: msa done in {:?} (elapsed {:?})",
                iter_idx + 1,
                iters,
                iter_start.elapsed(),
                start.elapsed()
            );
            thread::yield_now();
        }

        if let Some(exp) = expected.as_ref() {
            if &got != exp {
                return Err(abpoa::Error::InvalidInput(
                    format!(
                        "non-deterministic consensus within thread: expected={exp:?} got={got:?}"
                    )
                    .into(),
                ));
            }
        } else {
            expected = Some(got);
        }
    }

    if verbose {
        eprintln!("[{label}] done: total_elapsed={:?}", start.elapsed());
    }

    Ok(expected.unwrap_or_default())
}
