//! Aligns reads to a selected slice of the POA graph to demonstrate subgraph-restricted alignment

use abpoa::encode::encode_dna;
use abpoa::{
    Aligner, Alphabet, NodeId, OutputMode, Parameters, Result, SentinelNode, SubgraphRange,
};
use std::convert::TryFrom;

fn main() -> Result<()> {
    let sequences: Vec<&[u8]> = vec![
        b"CGTCAATCTATCGAAGCATACGCGGGCAGAGC",
        b"CGTCAATCTATCGAAGCATACGCGGCAGAGCC",
        b"CGTCAATCTATCGAAGCATACGCGGCAGAGTT",
    ];

    let mut aligner = Aligner::with_params(Parameters::new()?)?;

    let encoded: Vec<Vec<u8>> = sequences.iter().map(|seq| encode_dna(seq)).collect();
    let max_len = encoded.iter().map(Vec::len).max().unwrap_or(0);
    aligner.reset(max_len)?;
    let total_reads = i32::try_from(sequences.len()).unwrap();

    let whole_graph = SubgraphRange {
        beg: SentinelNode::Source.as_node_id(),
        end: SentinelNode::Sink.as_node_id(),
    };
    let first = aligner.align_sequence_to_subgraph(whole_graph, &encoded[0])?;
    aligner.add_subgraph_alignment(whole_graph, &encoded[0], &first, 0, total_reads, false)?;
    if first.cigar_len() > 0 {
        let graph = aligner.graph()?;
        let pretty = first.format_alignment(&graph, &encoded[0], Alphabet::Dna)?;
        println!("Seed read alignment:\n{}\n", pretty);
    }

    let include = SubgraphRange {
        beg: NodeId(2),
        end: NodeId((encoded[0].len() as i32 / 2) + 2),
    };
    let subgraph = aligner.subgraph_nodes(include.beg, include.end)?;
    for (idx, seq) in encoded.iter().enumerate().skip(1) {
        let alignment = aligner.align_sequence_to_subgraph(subgraph, seq)?;
        aligner.add_subgraph_alignment(
            subgraph,
            seq,
            &alignment,
            idx as i32,
            total_reads,
            false,
        )?;
        if alignment.cigar_len() > 0 {
            let graph = aligner.graph()?;
            let pretty = alignment.format_alignment(&graph, seq, Alphabet::Dna)?;
            println!("Alignment for read {}:\n{}\n", idx + 1, pretty);
        }
    }

    let result = aligner.finalize_msa(OutputMode::consensus_and_msa())?;
    println!("Consensus sequences:");
    for (idx, cluster) in result.clusters.iter().enumerate() {
        println!("  [{}] {}", idx, cluster.consensus);
    }

    println!("\nMSA rows:");
    for row in result.msa {
        println!("  {}", row);
    }

    Ok(())
}
