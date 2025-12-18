//! Demonstrates exporting an MSA and graph to disk, then restoring the graph from the saved files
//! (this may be a bit buggy if you use the default Parameters versus the custom Parameters used to save the graph).

use abpoa::{Aligner, Graph, OutputMode, Parameters, SequenceBatch};

fn main() -> abpoa::Result<()> {
    let seqs: Vec<&[u8]> = vec![
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"ATGCATGCATCGATCGATGCATCGATGCATGCATCGTTTTTTTTTTTTTTTTTTTATCGATCGATCG", // Reverse complement of long sequence
        b"CGATCGATCGATAAAAAAAAAAAAAAAAAAACGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"ATGCATGCATCGATCGATGCATCGATGCATGCATCGTTTTTTTTTTTTTTTTTTTATCGATCGATCG", // Reverse complement of long sequence
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
        b"ATGCATGCATCGATCGATGCATCGATGCATGCATCGATCGATCGATCG", // Reverse complement of short sequence
        b"CGATCGATCGATCGATGCATGCATCGATGCATCGATCGATGCATGCAT",
    ];

    let names: Vec<String> = seqs
        .iter()
        .enumerate()
        .map(|(i, _)| format!("read{}", i + 1))
        .collect();
    let name_refs: Vec<&str> = names.iter().map(|s| s.as_str()).collect();

    let mut params = Parameters::configure()?;
    params
        .set_max_consensus(2)?
        .set_use_read_ids(true)
        .set_ambiguous_strand(true); // Enable reverse complement handling
    let mut aligner = Aligner::with_params(params)?;

    let result = aligner.msa(SequenceBatch::from_sequences(&seqs).with_names(&name_refs))?;
    println!(
        "Consensus sequences: {:?}",
        result
            .clusters
            .iter()
            .map(|c| &c.consensus)
            .collect::<Vec<_>>()
    );

    let graph = aligner.graph()?;
    print_graph(&graph);

    // Write the MSA to a temporary FASTA file
    println!("\nWriting MSA to FASTA file, contents:");
    let mut buffer = Vec::new();
    aligner.write_msa_fasta(&mut buffer)?;
    println!("{}", String::from_utf8_lossy(&buffer));
    let temp_dir = std::env::temp_dir();
    let msa_path = temp_dir.join("sample_msa.fa");
    std::fs::write(&msa_path, &buffer)?;
    println!("Wrote MSA FASTA to: {}", msa_path.display());

    // Write the graph to a temporary GFA file
    println!("\nWriting graph to GFA file, contents:");
    let gfa_path = temp_dir.join("sample_graph.gfa");
    aligner.write_gfa_to_path(&gfa_path)?;
    let gfa_buffer = std::fs::read(&gfa_path)?;
    println!("{}", String::from_utf8_lossy(&gfa_buffer));
    println!("\nWrote graph to: {}", gfa_path.display());

    // Reload the graph from the saved FASTA or GFA file
    println!("\nReloading graph from FASTA or GFA file, contents:");
    // Note you can also use the Aligner::from_graph_file method to reload the graph from the MSA FASTA file but it will use default Parameters (n_consensus=1), probably good idea to add a Parameters argument to it
    // let mut reloaded_aligner = Aligner::from_graph_file(&msa_path, true)?;
    let mut params = Parameters::configure()?;
    params
        .set_max_consensus(2)?
        .set_incremental_graph_file(&gfa_path)?
        .set_use_read_ids(true)
        .set_outputs(OutputMode::CONSENSUS);
    let mut reloaded_aligner = Aligner::with_params(params)?;
    reloaded_aligner.restore_graph()?;

    let result = reloaded_aligner.finalize_msa()?;
    for (i, cluster) in result.clusters.iter().enumerate() {
        println!("Consensus {i}: {}", cluster.consensus);
    }

    std::fs::remove_file(&msa_path)?;
    std::fs::remove_file(&gfa_path)?;

    Ok(())
}

fn print_graph(graph: &Graph<'_>) {
    println!("Graph with {} sequences:", graph.sequences().len());
    for (idx, seq) in graph.sequences().iter().enumerate() {
        println!(
            "- #{idx}: {}, rc={} length={} bases",
            seq.name.as_str().unwrap(),
            seq.is_reverse_complement,
            seq.sequence.len()
        );
    }
}
