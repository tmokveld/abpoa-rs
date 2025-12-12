//! Demonstrates per-read subgraph ranges, aligning each sequence to different windows of the graph

use abpoa::encode::encode_dna;
use abpoa::{
    Aligner, Alphabet, ConsensusAlgorithm, NodeId, OutputMode, Parameters, Result, SentinelNode,
    SubgraphRange, Verbosity,
};

fn main() -> Result<()> {
    let first_n = 5; // Only consider the first n sequences (5 matches the example in the upstream C code)
    let seqs: Vec<&[u8]> = vec![
        b"CGTCAATCTATCGAAGCATACGCGGGCAGAGC",
        b"CCACGTCAATCTATCGAAGCATACGCGGCAGC",
        b"AATCTATCGAAGCATACG",
        b"CAATGCTAGTCGAAGCAGCTGCGGCAG",
        b"CAATGCTAGTCGAAGCAGCTGCGGCAG",
        b"CGTCAATCTATCGAAGCATTCTACGCGGCAGAGC",
        b"CGTCAATCTAGAAGCATACGCGGCAAGAGC",
        b"CGTCAATCTATCGGTAAAGCATACGCTCTGTAGC",
        b"CGTCAATCTATCTTCAAGCATACGCGGCAGAGC",
        b"CGTCAATGGATCGAGTACGCGGCAGAGC",
        b"CGTCAATCTAATCGAAGCATACGCGGCAGAGC",
    ];
    let beg_end_id: Vec<(i32, i32)> = vec![
        (0, 1),
        (2, 33),
        (6, 23),
        (5, 30),
        (5, 30),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
    ];

    let mut params = Parameters::configure()?;
    params
        .set_consensus(ConsensusAlgorithm::MostFrequent, 1, 0.0)?
        .set_inc_path_score(true)
        .set_sub_alignment(true)
        .set_verbosity(Verbosity::Info);

    let mut aligner = Aligner::with_params(params)?;
    let encoded: Vec<Vec<u8>> = seqs.iter().map(|s| encode_dna(s)).collect();
    let max_len = encoded
        .iter()
        .take(first_n)
        .map(Vec::len)
        .max()
        .unwrap_or(0);
    aligner.reset(max_len)?;
    let total_reads = seqs.iter().take(first_n).len() as i32;

    let whole = SubgraphRange {
        beg: SentinelNode::Source.as_node_id(),
        end: SentinelNode::Sink.as_node_id(),
    };

    for (idx, seq) in encoded.iter().take(first_n).enumerate() {
        let include = beg_end_id[idx];
        let range = if idx == 0 {
            println!("i: 0, beg: {}, end: {}", whole.beg.0, whole.end.0);
            whole
        } else {
            let range = aligner.subgraph_nodes(NodeId(include.0), NodeId(include.1))?;
            println!("i: {}, beg: {}, end: {}", idx, range.beg.0, range.end.0);
            range
        };
        let res = aligner.align_sequence_to_subgraph(range, seq)?;
        if res.cigar_len() > 0 {
            let graph = aligner.graph()?; // borrow immutably for pretty-printing
            let pretty = res.format_alignment(&graph, seq, Alphabet::Dna)?;
            println!("alignment for read {}:\n{}\n", idx + 1, pretty);
        }
        aligner.add_subgraph_alignment(range, seq, &res, idx as i32, total_reads, false)?;
    }

    let result = aligner.finalize_msa(OutputMode::consensus_and_msa())?;

    println!("\nMSA rows:");
    for row in result.msa {
        println!("  {}", row);
    }

    Ok(())
}
