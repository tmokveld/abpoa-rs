//! Generate two synthetic consensuses, simulate reads from both, mix them, and ask abPOA
//! to emit two consensus sequences

use abpoa::{Aligner, ConsensusAlgorithm, Parameters, Result, SequenceBatch};

#[path = "shared/sim.rs"]
mod sim;

fn main() -> Result<()> {
    let len_a = 200;
    let len_b = 100;

    let reads_per_consensus = 30;
    let seed_a = Some(42u64);
    let seed_b = Some(1337u64);

    let consensus_a = sim::make_consensus(sim::ConsensusOptions {
        len: len_a,
        seed: seed_a,
        ..Default::default()
    });

    // If None, generate B independently, and if Some(n), derive B from A by mutating every n-th
    // base and then resizing to len_b.
    // let mutate_every_b = None;
    let mutate_every_b = Some(5);

    let consensus_b = make_consensus_b(&consensus_a, len_b, seed_b, mutate_every_b);
    println!(
        "Generated consensus A length={}, consensus B length={}",
        consensus_a.len(),
        consensus_b.len()
    );

    let mut reads_a = sim::simulate_reads(
        &consensus_a,
        reads_per_consensus,
        sim::ErrorRates::new(0.002, 0.0025, 0.0025),
        Some(42),
    );
    let mut reads_b = sim::simulate_reads(
        &consensus_b,
        reads_per_consensus,
        sim::ErrorRates::new(0.002, 0.0025, 0.0025),
        Some(1337),
    );

    let mut reads: Vec<Vec<u8>> = Vec::with_capacity(reads_per_consensus * 2);
    for _ in 0..reads_per_consensus {
        reads.push(
            reads_a
                .next()
                .expect("reads_a should yield reads_per_consensus reads"),
        );
        reads.push(
            reads_b
                .next()
                .expect("reads_b should yield reads_per_consensus reads"),
        );
    }
    let read_refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();

    let mut params = Parameters::configure()?;
    params
        .set_consensus(ConsensusAlgorithm::HeaviestBundle, 2, 0.3)?
        .set_gap_on_right(true)
        .set_sub_alignment(true);
    let mut aligner = Aligner::with_params(params)?;

    let result = aligner.msa(SequenceBatch::from_sequences(&read_refs))?;

    result.clusters.iter().for_each(|cluster| {
        println!("Cluster read ids: {:?}", cluster.read_ids);
        println!("Cluster consensus: {:?}", cluster.consensus);
        println!("Cluster node ids: {:?}", cluster.node_ids);
        println!("Cluster coverage: {:?}", cluster.coverage);
        println!("Cluster phred: {:?}", cluster.phred);
    });

    println!(
        "Actual consensus A:\n{}",
        std::str::from_utf8(&consensus_a).unwrap()
    );
    println!(
        "Actual consensus B:\n{}",
        std::str::from_utf8(&consensus_b).unwrap()
    );

    println!("Inferred consensus clusters={}", result.clusters.len());
    for (idx, cluster) in result.clusters.iter().enumerate() {
        println!(
            "Inferred consensus {idx} length={}",
            cluster.consensus.len()
        );
        println!("Inferred consensus {idx}:\n{}", cluster.consensus);
    }
    assert_eq!(
        result.clusters.len(),
        2,
        "expected abPOA to emit two consensus sequences, got {}",
        result.clusters.len()
    );

    aligner.write_pog_to_path("test.pdf", &abpoa::PogDotOptions::default())?;

    Ok(())
}

fn mutate_consensus_every(consensus: &[u8], step: usize) -> Vec<u8> {
    assert!(step > 0);
    consensus
        .iter()
        .enumerate()
        .map(|(i, &b)| {
            if i % step != 0 {
                return b;
            }
            match b {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                b'T' => b'A',
                _ => b,
            }
        })
        .collect()
}

fn make_consensus_b(
    consensus_a: &[u8],
    len_b: usize,
    seed_b: Option<u64>,
    mutate_every: Option<usize>,
) -> Vec<u8> {
    match mutate_every {
        None => sim::make_consensus(sim::ConsensusOptions {
            len: len_b,
            seed: seed_b,
            ..Default::default()
        }),
        Some(mutate_every) => {
            let mutated = mutate_consensus_every(consensus_a, mutate_every);
            if len_b <= mutated.len() {
                mutated[..len_b].to_vec()
            } else {
                // Extend deterministically by cycling the mutated A-derived sequence.
                // (Keeps B "derived" from A while allowing arbitrary `len_b`.)
                let mut out = Vec::with_capacity(len_b);
                out.extend_from_slice(&mutated);
                for i in 0..(len_b - mutated.len()) {
                    out.push(mutated[i % mutated.len()]);
                }
                out
            }
        }
    }
}
