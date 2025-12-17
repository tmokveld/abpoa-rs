//! Show the difference between one-shot MSA (with and without seeding) and incremental MSA with minimizer seeding.
//! Only for the seeded one-shot run the minimizers are collected and the progressive guide tree is build and used.
//! The incremental C API does not use the minimizer seeding and progressive guide tree even if it is enabled in the parameters.

use abpoa::{Aligner, OutputMode, Parameters, Result, SequenceBatch, Verbosity};
use std::time::Instant;

fn main() -> Result<()> {
    let base = b"AAATTAGATGTCTACGAGTAAAGCGGTGAGTTTTACTAGCTTGAGGAATAAGAGCACAGAGTTCCAGGCTGACAGAAAAGAGGAACTGGGAAATTTGAATGACATGGGAGGACATCTCACACAACTGAAAGTCACAGAGAGAATATCACAGAGTAAAAATCTAAAATCAGCACTTCAACTTCATTCAAATATATGATGGCTGCTACATTTCACATTCATAAAAGGAGACTCCATAGAATCCAGCAGAAAACAACAGCTAAATTGCTAGTACAGAGCAGAGATTTCAACCATTGCATATAGCTCAGGAGTGAAAGTTTGGTGTTTGACTACGGAAAGAAGACTAATGTTAGAAAAGAGTCATTCTTCAAAGGAAAATAAACAAACCTATCTCTACAAAATA";
    let mut storage: Vec<Vec<u8>> = Vec::new();
    let options = [b'A', b'C', b'G', b'T'];
    for i in 0..30usize {
        let mut seq = base.repeat(1); // Low-complexity sequence kills the minimizer seeding (try setting repeat >1 to see the difference)
        // Mutate the sequence to have some variation across "reads".
        for pos in (i..seq.len()).step_by(30) {
            seq[pos] = options[(i + pos) % options.len()];
        }
        storage.push(seq);
    }
    let seqs: Vec<&[u8]> = storage.iter().map(|s| s.as_slice()).collect();

    println!("MSA one-shot without seeding - no minimizer logging should be expected!");
    let mut plain_params = Parameters::configure()?;
    plain_params
        .set_disable_seeding(true)
        .set_outputs(OutputMode::CONSENSUS)
        .set_verbosity(Verbosity::Info);
    let mut plain_aligner = Aligner::with_params(plain_params)?;

    let start = Instant::now();
    let result = plain_aligner.msa(SequenceBatch::from_sequences(&seqs))?;
    for (idx, cluster) in result.clusters.iter().enumerate() {
        println!("Consensus {}: {}", idx + 1, cluster.consensus);
    }

    println!("elapsed: {:?}\n", start.elapsed());

    println!(
        "MSA one-shot with minimizer seeding and progressive POA - minimizer logging should be expected!"
    );
    let mut seeded_params = Parameters::configure()?;
    seeded_params
        .set_minimizer_seeding(6, 4, 5)?
        .set_progressive_poa(true)
        .set_outputs(OutputMode::CONSENSUS)
        .set_verbosity(Verbosity::Info);
    let mut seeded_aligner = Aligner::with_params(seeded_params)?;

    let start = Instant::now();
    let result = seeded_aligner.msa(SequenceBatch::from_sequences(&seqs))?;
    for (idx, cluster) in result.clusters.iter().enumerate() {
        println!("Consensus {}: {}", idx + 1, cluster.consensus);
    }
    println!("elapsed: {:?}\n", start.elapsed());

    println!(
        "Incremental build with the same parameters - no minimizer or guide-tree logging should be expected since this code-path is not used by the incremental API"
    );
    let mut inc_params = Parameters::configure()?;
    inc_params
        .set_minimizer_seeding(7, 4, 10)?
        .set_progressive_poa(true)
        .set_outputs(OutputMode::CONSENSUS)
        .set_verbosity(Verbosity::Info);
    let mut incremental_aligner = Aligner::with_params(inc_params)?;

    let start = Instant::now();
    incremental_aligner.msa_in_place(SequenceBatch::from_sequences(&seqs))?;
    let result = incremental_aligner.finalize_msa()?;
    for (idx, cluster) in result.clusters.iter().enumerate() {
        println!("Consensus {}: {}", idx + 1, cluster.consensus);
    }
    println!("elapsed: {:?}\n", start.elapsed());

    Ok(())
}
