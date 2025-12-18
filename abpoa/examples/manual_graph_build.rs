//! Compares the high-level msa helper with an equivalent manual graph-building flow.

use abpoa::{
    AlignMode, Aligner, Alphabet, ConsensusAlgorithm, Parameters, Result, Scoring, SequenceBatch,
    Verbosity, encode::encode_dna,
};

fn main() -> Result<()> {
    let sequences: Vec<&[u8]> = vec![
        b"ACGTGTACAGTTGAC",
        b"AGGTACACGTTAC",
        b"AGTGTCACGTTGAC",
        b"ACGTGTACATTGAC",
    ];

    let truth = vec![
        "ACGTGTACA-GTTGAC",
        "A-G-GTACACGTT-AC",
        "A-GTGT-CACGTTGAC",
        "ACGTGTACA--TTGAC",
    ];
    let truth_consensus = "ACGTGTACACGTTGAC";

    let mut params = Parameters::configure()?;
    params
        .set_align_mode(AlignMode::Global)
        .set_scoring_scheme(Scoring::convex(2, 4, 4, 2, 24, 1))?
        .set_consensus(ConsensusAlgorithm::HeaviestBundle, 1, 0.0)?
        .set_verbosity(Verbosity::Info);

    let mut aligner = Aligner::with_params(params)?;
    let result = aligner.msa(SequenceBatch::from_sequences(&sequences))?;

    assert_eq!(result.msa.len(), truth.len());
    for (row, expected) in result.msa.iter().zip(truth.iter()) {
        assert_eq!(row, expected);
    }

    let consensus = result
        .clusters
        .first()
        .map(|c| c.consensus.as_str())
        .unwrap();
    assert_eq!(consensus, truth_consensus);

    // The equivalent manual graph-building flow using the low-level API:
    let mut manual_params = Parameters::configure()?;
    manual_params
        .set_align_mode(AlignMode::Global)
        .set_scoring_scheme(Scoring::convex(2, 4, 4, 2, 24, 1))?
        .set_consensus(ConsensusAlgorithm::HeaviestBundle, 1, 0.0)?
        .set_verbosity(Verbosity::None);
    let mut manual = Aligner::with_params(manual_params)?;
    let encoded: Vec<Vec<u8>> = sequences.iter().map(|s| encode_dna(s)).collect();
    let max_len = encoded.iter().map(Vec::len).max().unwrap_or(0);
    manual.reset(max_len)?; // Reset the aligner to the maximum length of the sequences
    let total_reads = sequences.len() as i32;
    for (idx, seq) in encoded.iter().enumerate() {
        let aln = manual.align_sequence_raw(seq)?; // Align the sequence
        if aln.cigar_len() > 0 {
            let graph = manual.graph()?; // Render per-read alignment before mutating the graph
            let pretty = aln.format_alignment(&graph, seq, Alphabet::Dna)?;
            println!("Read {} alignment:\n{}\n", idx + 1, pretty);
        }
        manual.add_alignment(seq, &aln, idx as i32, total_reads)?; // Add the alignment to the graph
    }
    let manual_result = manual.finalize_msa()?; // Finalize the MSA

    assert_eq!(manual_result.msa.len(), truth.len());
    for (row, expected) in manual_result.msa.iter().zip(truth.iter()) {
        assert_eq!(row, expected);
    }

    let manual_consensus = manual_result
        .clusters
        .first()
        .map(|c| c.consensus.as_str())
        .unwrap();
    assert_eq!(manual_consensus, truth_consensus);

    Ok(())
}
