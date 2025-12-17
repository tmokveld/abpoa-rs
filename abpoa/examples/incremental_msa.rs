//! Aligns sequences in multiple passes to show how the POA graph grows before generating a final MSA.

use abpoa::{Aligner, Parameters, Result, SequenceBatch};

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
        b"CGATCGATCGATGGGGGGGGGGGGGGGGGGGCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
    ];
    let first_part = 5;
    let second_part = 3;
    let third_part = seqs.len() - first_part - second_part;

    let mut params = Parameters::configure()?;
    params.set_max_consensus(2)?;
    let mut aligner = Aligner::with_params(params)?;

    // Align the first few sequences
    aligner.msa_in_place(SequenceBatch::from_sequences(&seqs[..first_part]))?;
    let initial = aligner.finalize_msa()?;
    println!("Initial MSA (first {first_part} sequences):");
    println!("\nMSA rows:");
    for row in &initial.msa {
        println!("  {}", row);
    }
    let graph = aligner.graph()?;
    println!(
        "\nExisting graph has {} sequences, {} nodes\n",
        initial.msa.len(),
        graph.node_count(),
    );

    aligner.write_pog_to_path("test_1_initial.pdf", &abpoa::PogDotOptions::default())?;

    println!(
        "Appending the second part of {} sequences one by one...",
        second_part
    );
    // Append sequences one at a time without generating output
    for seq in seqs[first_part..first_part + second_part].iter() {
        aligner.add_sequences(SequenceBatch::from_sequences(&[seq]))?;
    }

    let secondary = aligner.finalize_msa()?;
    println!("Secondary MSA (after appending the second part of {second_part} sequences):");
    println!("\nMSA rows:");
    for row in &secondary.msa {
        println!("  {}", row);
    }

    let graph = aligner.graph()?;
    println!(
        "After appending the second part, existing graph has {} sequences, {} nodes\n",
        first_part + second_part,
        graph.node_count(),
    );

    aligner.write_pog_to_path("test_2_second_part.pdf", &abpoa::PogDotOptions::default())?;

    println!(
        "Appending the remaining {} sequences together...",
        third_part
    );
    // Append the remaining sequences together and output the MSA and consensus
    if third_part > 0 {
        aligner.add_sequences(SequenceBatch::from_sequences(
            &seqs[seqs.len() - third_part..],
        ))?;
    }

    let result = aligner.finalize_msa()?;
    let graph = aligner.graph()?;
    println!(
        "> Final MSA after all sequences ({} total, {} nodes):",
        result.msa.len(),
        graph.node_count()
    );

    aligner.write_pog_to_path("test_3_final.pdf", &abpoa::PogDotOptions::default())?;

    println!("\nMSA rows:");
    for row in result.msa {
        println!("  {}", row);
    }

    Ok(())
}
