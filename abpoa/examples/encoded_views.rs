//! Demonstrate zero-copy encoded output views (`EncodedMsaView` / `EncodedMsaRows` / `EncodedClusterView`).

use abpoa::{Aligner, EncodedClusterView, EncodedMsaRows, Parameters, SequenceBatch, encode};

fn main() -> abpoa::Result<()> {
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

    let mut params = Parameters::configure()?;
    params.set_max_consensus(2)?;
    let mut aligner = Aligner::with_params(params)?;
    let view = aligner.msa_view_encoded(SequenceBatch::from_sequences(&seqs))?;

    println!(
        "alphabet: {:?}\nsequences: {}\nmsa_len: {}\nclusters: {}",
        view.alphabet(),
        view.sequence_count(),
        view.msa_len(),
        view.cluster_count()
    );

    // Random-access to a single MSA row (still borrowed, zero-copy)
    if let Some(row0) = view.msa_row(1) {
        println!("\nmsa_row(1): {}", encode::decode_dna(row0));
    }

    // Iterate MSA rows with `EncodedMsaRows`
    let mut rows: EncodedMsaRows<'_> = view.msa_rows();
    println!("\nFirst 3 MSA rows (decoded):");
    for (idx, row) in rows.by_ref().take(3).enumerate() {
        println!("  row {idx}: {}", encode::decode_dna(row));
    }
    println!("  ... ({} more rows)", rows.len());

    // Random-access to a single cluster view
    if let Some(cluster0) = view.cluster(1) {
        print_cluster("cluster(1)", cluster0);
    }

    // Iterate consensus clusters with `EncodedClusterView`
    println!("\nAll clusters:");
    for (idx, cluster) in view.clusters().enumerate() {
        print_cluster(&format!("cluster {idx}"), cluster);
    }

    Ok(())
}

fn print_cluster(label: &str, cluster: EncodedClusterView<'_>) {
    println!("\n{label}:");
    println!("  read ids: {:?}", cluster.read_ids());
    println!("  consensus: {}", encode::decode_dna(cluster.consensus()));
    println!("  node_ids: {:?}", &cluster.node_ids_raw());
}
