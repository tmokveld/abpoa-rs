//! Rust port of the upstream `example.c`, showing quality-weighted MSA and graph/consensus outputs

use abpoa::{Aligner, Parameters, Result, SequenceBatch, Verbosity};

fn main() -> Result<()> {
    let seqs: Vec<&[u8]> = vec![
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
    ];
    let weights = quality_weights(&seqs);
    let weight_refs: Vec<&[i32]> = weights.iter().map(|w| w.as_slice()).collect();

    let mut params = Parameters::configure()?;
    params
        .set_minimizer_seeding(9, 6, 10)?
        .set_progressive_poa(false)
        .set_max_consensus(2)?
        .set_use_quality(true)
        .set_verbosity(Verbosity::Info);

    let mut aligner = Aligner::with_params(params)?;
    let result =
        aligner.msa(SequenceBatch::from_sequences(&seqs).with_quality_weights(&weight_refs))?;

    println!("Consensus sequences:");
    for cluster in result.clusters {
        println!("Consensus: {}", cluster.consensus);
    }

    println!("=== msa using manual print ===");
    println!("\nMSA rows:");
    for row in result.msa {
        println!("  {}", row);
    }

    println!("\n=== msa using write_msa_fasta ===");
    let mut stdout = std::io::stdout();
    aligner.write_msa_fasta(&mut stdout)?;

    println!("\n=== gfa using write_gfa ===");
    aligner.write_gfa(&mut stdout)?;

    println!("\n=== consensus using write_consensus_fasta ===");
    aligner.write_consensus_fasta(&mut stdout)?;

    aligner.write_pog_to_path("test.pdf", abpoa::PogDotOptions::default())?;

    Ok(())
}

// Mirror the quality weights function in the example.c file
fn quality_weights(seqs: &[&[u8]]) -> Vec<Vec<i32>> {
    seqs.iter()
        .map(|seq| {
            seq.iter()
                .enumerate()
                .map(|(idx, _)| if idx >= 12 { 2 } else { 0 })
                .collect()
        })
        .collect()
}
