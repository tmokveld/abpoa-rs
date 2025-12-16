use super::*;
use crate::encode::{decode_aa, decode_dna, encode_dna};
use crate::result::EncodedCluster;
use crate::{EncodedMsaResult, EncodedMsaView, MsaResult, OutputMode, SentinelNode};
use std::{
    process,
    time::{SystemTime, UNIX_EPOCH},
};

#[test]
fn aligner_init_and_free() {
    let mut aligner = Aligner::new().unwrap();
    assert!(!aligner.as_mut_ptr().is_null());
}

#[test]
fn from_graph_file_restores_graph_and_toggle_read_ids() {
    let mut aligner = Aligner::new().unwrap();
    aligner
        .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
        .unwrap();

    let path = std::env::temp_dir().join(format!(
        "abpoa_restore_{}_{}.gfa",
        process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    aligner.write_gfa_to_path(&path).unwrap();

    let mut restored = Aligner::from_graph_file(&path, false).unwrap();
    let graph = restored.graph().unwrap();
    assert!(
        !graph.is_empty(),
        "restored graph should include aligned nodes"
    );
    let raw_params = unsafe { restored.params_mut().as_mut_ptr().as_ref() }.unwrap();
    assert!(
        raw_params.use_read_ids() == 0,
        "preserve_read_ids flag should propagate to parameters"
    );

    std::fs::remove_file(&path).unwrap_or_default();
}

#[test]
fn raw_alignment_exposes_cigar_and_counts() {
    let mut aligner = Aligner::new().unwrap();
    {
        let params = aligner.params_mut();
        params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
        params.set_alphabet(Alphabet::Dna).unwrap();
    }

    let ref_seq = encode_dna(b"ACGT");
    let query = encode_dna(b"AGGT");
    let max_len = ref_seq.len().max(query.len());
    aligner.reset(max_len).unwrap();

    let first = aligner.align_sequence_raw(&ref_seq).unwrap();
    aligner.add_alignment(&ref_seq, &first, 0, 2).unwrap();

    let second = aligner.align_sequence_raw(&query).unwrap();

    assert_eq!(second.aligned_bases(), query.len() as i32);
    assert!(second.matched_bases() < second.aligned_bases());

    let ops: Vec<_> = second.cigar().collect();
    assert_eq!(ops.len() as i32, second.cigar_len());
    assert!(
        ops.iter()
            .any(|op| matches!(op, GraphCigarOp::Aligned { .. })),
        "graph cigar should contain aligned nodes"
    );
}

#[test]
fn dna_msa_variants() {
    let sequences = [b"ACGT".as_ref(), b"ACGT".as_ref(), b"ACGG".as_ref()];
    let mut aligner = Aligner::new().unwrap();
    let decoded = aligner
        .msa(
            SequenceBatch::from_sequences(&sequences),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    assert_eq!(decoded.msa.len(), sequences.len());
    assert_eq!(decoded.msa[0].len(), 4);
    assert!(!decoded.clusters.is_empty());
    assert_eq!(decoded.clusters[0].consensus, "ACGT");
    for cluster in &decoded.clusters {
        let len = cluster.consensus.len();
        assert_eq!(cluster.node_ids.len(), len);
        assert_eq!(cluster.coverage.len(), len);
        assert_eq!(cluster.phred.len(), len);
        assert!(
            cluster.phred.iter().all(|q| q.is_ascii_graphic()),
            "phred scores should be printable FASTQ bytes"
        );
    }

    let encoded = aligner
        .msa_encoded(
            SequenceBatch::from_sequences(&sequences),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    assert_eq!(encoded.msa.len(), decoded.msa.len());
    for (enc_row, dec_row) in encoded.msa.iter().zip(decoded.msa.iter()) {
        assert_eq!(decode_dna(enc_row), *dec_row);
    }

    assert_eq!(encoded.clusters.len(), decoded.clusters.len());
    for (enc_cluster, dec_cluster) in encoded.clusters.iter().zip(decoded.clusters.iter()) {
        assert_eq!(enc_cluster.read_ids, dec_cluster.read_ids);
        assert_eq!(decode_dna(&enc_cluster.consensus), dec_cluster.consensus);
        assert_eq!(enc_cluster.node_ids, dec_cluster.node_ids);
        assert_eq!(enc_cluster.coverage, dec_cluster.coverage);
        assert_eq!(enc_cluster.phred, dec_cluster.phred);
    }

    let sequences_with_gap = [b"ACGT".as_ref(), b"AGT".as_ref()];
    let mut aligner = Aligner::new().unwrap();
    let gap_result = aligner
        .msa(
            SequenceBatch::from_sequences(&sequences_with_gap),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    assert_eq!(gap_result.msa.len(), sequences_with_gap.len());
    assert!(
        gap_result
            .msa
            .iter()
            .all(|row| row.len() == gap_result.msa[0].len()),
        "msa rows should be padded to the same length"
    );
    assert!(
        gap_result.msa.iter().any(|row| row.contains('-')),
        "msa should surface gap sentinels as '-'"
    );

    let sequences_quality = [b"ACGT".as_ref(), b"ACGT".as_ref()];
    let qual_a = [30, 30, 30, 30];
    let qual_b = [10, 10, 10, 10];
    let qualities: [&[i32]; 2] = [qual_a.as_slice(), qual_b.as_slice()];
    let mut aligner = Aligner::new().unwrap();

    let qual_result = aligner
        .msa(
            SequenceBatch::from_sequences(&sequences_quality).with_quality_weights(&qualities),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    assert_eq!(qual_result.msa.len(), sequences_quality.len());
    assert!(!qual_result.clusters.is_empty());
}

#[test]
fn consensus_writer_matches_alignment() {
    let mut aligner = Aligner::new().unwrap();
    aligner
        .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGT"]))
        .unwrap();

    let mut buffer = Vec::new();
    aligner.write_consensus_fasta(&mut buffer).unwrap();
    let output = String::from_utf8(buffer).unwrap();
    assert!(
        output.contains(">Consensus_sequence"),
        "consensus header should be present"
    );
    assert!(
        output.contains("ACGT"),
        "consensus sequence should match aligned reads"
    );
}

#[test]
fn consensus_fastq_writer_produces_valid_records() {
    let mut aligner = Aligner::new().unwrap();
    aligner
        .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGT"]))
        .unwrap();

    let mut buffer = Vec::new();
    aligner.write_consensus_fastq(&mut buffer).unwrap();
    let output = String::from_utf8(buffer).unwrap();
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        lines.len() % 4,
        0,
        "FASTQ output should have 4-line records"
    );
    for record in lines.chunks(4) {
        assert!(
            record[0].starts_with("@Consensus_sequence"),
            "FASTQ record should have @ header"
        );
        assert!(
            record[2].starts_with("+Consensus_sequence"),
            "FASTQ record should have + header"
        );
        assert_eq!(
            record[1].len(),
            record[3].len(),
            "FASTQ record should have matching sequence and quality lengths"
        );
        assert!(
            record[3].as_bytes().iter().all(|q| q.is_ascii_graphic()),
            "FASTQ quality line should contain printable bytes"
        );
    }
}

#[test]
fn msa_writer_includes_names_and_reverse_complements() {
    let mut params = Parameters::configure().unwrap();
    params.set_ambiguous_strand(true);
    let mut aligner = Aligner::with_params(params).unwrap();
    let sequences = [b"ACGTTGC".as_ref(), b"GCAACGT".as_ref()]; // read2 is the RC of read1
    let names = ["read1", "read2"];

    aligner
        .msa(
            SequenceBatch::from_sequences(&sequences).with_names(&names),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    let mut buffer = Vec::new();
    aligner.write_msa_fasta(&mut buffer).unwrap();
    let output = String::from_utf8(buffer).unwrap();
    assert!(output.contains(">read1"), "msa should include read1 header");
    assert!(output.contains(">read2"), "msa should include read2 header");
    assert!(
        output.contains(">Consensus_sequence"),
        "msa writer should include consensus rows"
    );
    assert!(
        output.contains("_reverse_complement"),
        "msa output should flag reverse complements when ambiguous strand is enabled"
    );
}

#[test]
fn gfa_writer_generates_header() {
    let mut aligner = Aligner::new().unwrap();
    aligner
        .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
        .unwrap();

    let mut buffer = Vec::new();
    aligner.write_gfa(&mut buffer).unwrap();
    let output = String::from_utf8(buffer).unwrap();
    assert!(
        output.contains("H\tVN:Z:1.0"),
        "GFA header should be present"
    );
    assert!(
        output.contains("\nS\t"),
        "GFA output should include at least one segment"
    );
}

#[test]
fn graph_alignment_round_trip() {
    let mut aligner = Aligner::new().unwrap();
    {
        let params = aligner.params_mut();
        params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
        params.set_alphabet(Alphabet::Dna).unwrap();
    }

    let seqs = [encode_dna(b"ACGT"), encode_dna(b"ACGG")];
    let max_len = seqs.iter().map(Vec::len).max().unwrap();
    aligner.reset(max_len).unwrap();

    let first = aligner.align_sequence_raw(&seqs[0]).unwrap();
    aligner.add_alignment(&seqs[0], &first, 0, 2).unwrap();

    let second = aligner.align_sequence_raw(&seqs[1]).unwrap();
    aligner.add_alignment(&seqs[1], &second, 1, 2).unwrap();

    let params_ptr = {
        let params = aligner.params_mut();
        params.as_mut_ptr()
    };
    unsafe { sys::abpoa_generate_consensus(aligner.as_mut_ptr(), params_ptr) };
    let abc = unsafe { (*aligner.as_ptr()).abc };
    let result = unsafe { MsaResult::from_raw(abc, Alphabet::Dna) };
    unsafe { sys::abpoa_clean_msa_cons(aligner.as_mut_ptr()) };

    let mut expected_aligner = Aligner::new().unwrap();
    let expected = expected_aligner
        .msa(
            SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    assert_eq!(result.clusters.len(), 1);
    assert_eq!(expected.clusters.len(), 1);
    assert_eq!(
        result.clusters[0].consensus, expected.clusters[0].consensus,
        "incremental consensus should match one-shot run"
    );
}

#[test]
fn incremental_api_matches_one_shot() {
    let mut aligner = Aligner::new().unwrap();
    aligner
        .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
        .unwrap();
    let incremental = aligner
        .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
        .unwrap();

    let mut expected_aligner = Aligner::new().unwrap();
    let expected = expected_aligner
        .msa(
            SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    assert_eq!(incremental.clusters.len(), expected.clusters.len());
    assert_eq!(
        incremental.clusters[0].consensus, expected.clusters[0].consensus,
        "incremental consensus should match one-shot run"
    );
    assert_eq!(incremental.msa.len(), expected.msa.len());
}

#[test]
fn finalize_msa_recomputes_after_cleanup() {
    let mut aligner = Aligner::new().unwrap();
    aligner
        .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGT"]))
        .unwrap();

    let first = aligner
        .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
        .unwrap();
    assert!(!first.clusters.is_empty(), "initial consensus should exist");

    let second = aligner
        .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
        .unwrap();
    assert_eq!(second.clusters.len(), first.clusters.len());
    assert_eq!(second.clusters[0].consensus, first.clusters[0].consensus);
    assert_eq!(second.msa.len(), first.msa.len());
}

#[test]
fn adding_sequences_updates_graph() {
    let mut aligner = Aligner::new().unwrap();
    aligner
        .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
        .unwrap();
    aligner
        .add_sequences(SequenceBatch::from_sequences(&[b"ACGA"]))
        .unwrap();

    let incremental = aligner
        .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
        .unwrap();

    let mut expected_aligner = Aligner::new().unwrap();
    let expected = expected_aligner
        .msa(
            SequenceBatch::from_sequences(&[b"ACGT", b"ACGG", b"ACGA"]),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    assert_eq!(incremental.msa.len(), expected.msa.len());
    assert_eq!(
        incremental.clusters[0].consensus, expected.clusters[0].consensus,
        "incremental consensus should match one-shot run"
    );
}

#[test]
fn subgraph_alignment_matches_one_shot() {
    let mut aligner = Aligner::new().unwrap();
    {
        let params = aligner.params_mut();
        params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
        params.set_alphabet(Alphabet::Dna).unwrap();
    }

    let seqs = [
        encode_dna(b"ACGT"),
        encode_dna(b"ACGG"),
        encode_dna(b"ACGA"),
    ];
    let max_len = seqs.iter().map(Vec::len).max().unwrap();
    aligner.reset(max_len).unwrap();

    let total_reads = to_i32(seqs.len(), "too many sequences for abpoa").unwrap();
    let whole_graph = SubgraphRange {
        beg: SentinelNode::Source.as_node_id(),
        end: SentinelNode::Sink.as_node_id(),
    };

    let first = aligner
        .align_sequence_to_subgraph(whole_graph, &seqs[0])
        .unwrap();
    aligner
        .add_subgraph_alignment(whole_graph, &seqs[0], &first, 0, total_reads, false)
        .unwrap();

    let last_node = to_i32(seqs[0].len() + 1, "sequence length exceeds i32").unwrap();
    let range = aligner
        .subgraph_nodes(NodeId(2), NodeId(last_node))
        .unwrap();
    assert_eq!(range.beg, SentinelNode::Source.as_node_id());
    assert_eq!(range.end, SentinelNode::Sink.as_node_id());

    let second = aligner.align_sequence_to_subgraph(range, &seqs[1]).unwrap();
    aligner
        .add_subgraph_alignment(range, &seqs[1], &second, 1, total_reads, false)
        .unwrap();

    let third = aligner.align_sequence_to_subgraph(range, &seqs[2]).unwrap();
    aligner
        .add_subgraph_alignment(range, &seqs[2], &third, 2, total_reads, false)
        .unwrap();

    let incremental = aligner
        .finalize_msa(OutputMode::CONSENSUS | OutputMode::MSA)
        .unwrap();

    let mut expected_aligner = Aligner::new().unwrap();
    let expected = expected_aligner
        .msa(
            SequenceBatch::from_sequences(&[b"ACGT", b"ACGG", b"ACGA"]),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    assert_eq!(
        incremental.clusters[0].consensus,
        expected.clusters[0].consensus
    );
    assert_eq!(incremental.msa.len(), expected.msa.len());
    assert_eq!(incremental.msa[0].len(), expected.msa[0].len());
}

#[test]
fn format_alignment_renders_indels() {
    fn check_case<F>(reference: &[u8], query: &[u8], expected_pretty: &str, expects_op: F)
    where
        F: Fn(&GraphCigarOp) -> bool,
    {
        let mut aligner = Aligner::new().unwrap();
        {
            let params = aligner.params_mut();
            params.set_outputs(OutputMode::CONSENSUS | OutputMode::MSA);
            params.set_alphabet(Alphabet::Dna).unwrap();
        }

        let reference = encode_dna(reference);
        let query = encode_dna(query);
        let max_len = reference.len().max(query.len());
        aligner.reset(max_len).unwrap();

        let base_alignment = aligner.align_sequence_raw(&reference).unwrap();
        aligner
            .add_alignment(&reference, &base_alignment, 0, 2)
            .unwrap();

        let query_alignment = aligner.align_sequence_raw(&query).unwrap();
        let ops: Vec<_> = query_alignment.cigar().collect();
        assert!(
            ops.iter().any(expects_op),
            "graph cigar should include expected indel"
        );

        let graph = aligner.graph().unwrap();
        let pretty = query_alignment
            .format_alignment(&graph, &query, Alphabet::Dna)
            .unwrap();
        assert_eq!(pretty, expected_pretty);
    }

    check_case(b"ACGT", b"AACGT", "-ACGT\n ||||\nAACGT", |op| {
        matches!(
            op,
            GraphCigarOp::QueryRun {
                op: CigarOp::Insertion,
                len,
                ..
            } if *len >= 1
        )
    });

    check_case(
        b"ACGT",
        b"ACG",
        "ACGT\n||| \nACG-",
        |op| matches!(op, GraphCigarOp::Deletion { len, .. } if *len >= 1),
    );
}

#[test]
fn msa_encoded_matches_decoded_for_amino_acids() {
    let mut aligner = Aligner::new().unwrap();
    {
        let params = aligner.params_mut();
        params.set_alphabet(Alphabet::AminoAcid).unwrap();
    }

    let sequences = [b"ACDE".as_ref(), b"ACDF".as_ref()];

    let decoded = aligner
        .msa(
            SequenceBatch::from_sequences(&sequences),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();
    assert_eq!(decoded.msa.len(), sequences.len());
    assert!(
        !decoded.clusters.is_empty(),
        "consensus cluster should be returned for amino acid input"
    );
    assert!(
        decoded.clusters[0].consensus.starts_with("ACD"),
        "consensus should decode amino acid alphabet"
    );

    let encoded = aligner
        .msa_encoded(
            SequenceBatch::from_sequences(&sequences),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();

    assert_eq!(encoded.msa.len(), decoded.msa.len());
    for (enc_row, dec_row) in encoded.msa.iter().zip(decoded.msa.iter()) {
        assert_eq!(decode_aa(enc_row), *dec_row);
    }

    assert_eq!(encoded.clusters.len(), decoded.clusters.len());
    for (enc_cluster, dec_cluster) in encoded.clusters.iter().zip(decoded.clusters.iter()) {
        assert_eq!(enc_cluster.read_ids, dec_cluster.read_ids);
        assert_eq!(decode_aa(&enc_cluster.consensus), dec_cluster.consensus);
        assert_eq!(enc_cluster.node_ids, dec_cluster.node_ids);
        assert_eq!(enc_cluster.coverage, dec_cluster.coverage);
        assert_eq!(enc_cluster.phred, dec_cluster.phred);
    }
}

fn copy_encoded_view(view: &EncodedMsaView<'_>) -> EncodedMsaResult {
    let msa = view.msa_rows().map(|row| row.to_vec()).collect();
    let clusters = view
        .clusters()
        .map(|cluster| EncodedCluster {
            read_ids: cluster.read_ids().to_vec(),
            consensus: cluster.consensus().to_vec(),
            node_ids: cluster.node_ids_raw().iter().copied().map(NodeId).collect(),
            coverage: cluster.coverage().to_vec(),
            phred: cluster
                .phred_scores_raw()
                .iter()
                .copied()
                .map(|score| u8::try_from(score).unwrap_or(b'!'))
                .collect(),
        })
        .collect();
    EncodedMsaResult { msa, clusters }
}

#[test]
fn msa_view_encoded_matches_owned_one_shot() {
    let sequences = [b"ACGT".as_ref(), b"ACGG".as_ref(), b"ACG".as_ref()];
    let mut aligner = Aligner::new().unwrap();

    let from_view = {
        let view = aligner
            .msa_view_encoded(
                SequenceBatch::from_sequences(&sequences),
                OutputMode::CONSENSUS | OutputMode::MSA,
            )
            .unwrap();
        assert_eq!(view.sequence_count(), sequences.len());
        assert!(view.msa_len() > 0);
        assert!(view.consensus_msa_row(0).is_some());
        copy_encoded_view(&view)
    };

    let owned = aligner
        .msa_encoded(
            SequenceBatch::from_sequences(&sequences),
            OutputMode::CONSENSUS | OutputMode::MSA,
        )
        .unwrap();
    assert_eq!(from_view, owned);
}

#[test]
fn finalize_msa_view_encoded_matches_owned() {
    let sequences = [b"ACGT".as_ref(), b"ACGG".as_ref(), b"ACG".as_ref()];
    let mut aligner = Aligner::new().unwrap();
    aligner
        .msa_in_place(SequenceBatch::from_sequences(&sequences))
        .unwrap();

    let from_view = {
        let view = aligner
            .finalize_msa_view_encoded(OutputMode::CONSENSUS | OutputMode::MSA)
            .unwrap();
        assert_eq!(view.sequence_count(), sequences.len());
        assert!(view.msa_len() > 0);
        assert!(view.consensus_msa_row(0).is_some());
        copy_encoded_view(&view)
    };

    let owned = aligner
        .finalize_msa_encoded(OutputMode::CONSENSUS | OutputMode::MSA)
        .unwrap();
    assert_eq!(from_view, owned);
}

#[test]
fn read_id_less_graph_rejects_msa_gfa_and_multi_consensus() {
    let mut aligner = Aligner::new().unwrap();
    {
        let params = aligner.params_mut();
        params.set_outputs(OutputMode::CONSENSUS);
        params.set_use_read_ids(false);
    }

    aligner
        .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
        .unwrap();

    // Single-consensus is still supported without read ids.
    let consensus = aligner.finalize_msa(OutputMode::CONSENSUS).unwrap();
    assert!(!consensus.clusters.is_empty());

    assert!(matches!(
        aligner.finalize_msa(OutputMode::MSA),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        aligner.write_msa_fasta(&mut Vec::new()),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        aligner.write_gfa(&mut Vec::new()),
        Err(Error::InvalidInput(_))
    ));

    {
        let params = aligner.params_mut();
        params.set_max_consensus(2).unwrap();
    }
    assert!(matches!(
        aligner.finalize_msa(OutputMode::CONSENSUS),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn can_generate_msa_and_multi_consensus_after_consensus_only_build_with_read_ids() {
    let mut params = Parameters::configure().unwrap();
    params.set_outputs(OutputMode::CONSENSUS);
    params.set_use_read_ids(true);

    let mut aligner = Aligner::with_params(params).unwrap();
    aligner
        .msa_in_place(SequenceBatch::from_sequences(&[b"ACGT", b"ACGG"]))
        .unwrap();

    let mut msa = Vec::new();
    aligner.write_msa_fasta(&mut msa).unwrap();
    assert!(!msa.is_empty());

    aligner.params_mut().set_max_consensus(2).unwrap();
    let mut consensus = Vec::new();
    aligner.write_consensus_fasta(&mut consensus).unwrap();
    assert!(!consensus.is_empty());
}

#[test]
fn manual_graph_editing_and_topology_sorting() {
    let mut aligner = Aligner::new().unwrap();
    {
        let params = aligner.params_mut();
        params.set_alphabet(Alphabet::Dna).unwrap();
    }

    let base = encode_dna(b"A")[0];
    let src = SentinelNode::Source.as_node_id();
    let sink = SentinelNode::Sink.as_node_id();
    let mid = aligner.add_node(base).unwrap();
    aligner.add_edge(src, mid, 1, false).unwrap();
    aligner.add_edge(mid, sink, 1, false).unwrap();

    aligner.ensure_topological().unwrap();
    let graph = aligner.graph().unwrap();
    assert_eq!(graph.node_count(), 3);

    let src_idx = graph.node_index(src).unwrap();
    let mid_idx = graph.node_index(mid).unwrap();
    let sink_idx = graph.node_index(sink).unwrap();
    assert!(src_idx < mid_idx && mid_idx < sink_idx);

    let weight = graph.node_weight(src).unwrap();
    assert_eq!(weight, 1);
}

#[test]
fn refresh_helpers_require_allocated_buffers() {
    let mut aligner = Aligner::new().unwrap();
    let src = SentinelNode::Source.as_node_id();
    let sink = SentinelNode::Sink.as_node_id();
    let mid = aligner.add_node(encode_dna(b"C")[0]).unwrap();
    aligner.add_edge(src, mid, 1, false).unwrap();
    aligner.add_edge(mid, sink, 1, false).unwrap();

    assert!(
        matches!(
            aligner.refresh_node_indices(),
            Err(Error::NullPointer(_) | Error::InvalidInput(_))
        ),
        "refresh should fail until topology buffers exist"
    );
    assert!(
        matches!(
            aligner.refresh_node_remaining(),
            Err(Error::NullPointer(_) | Error::InvalidInput(_))
        ),
        "refresh should fail until topology buffers exist"
    );

    aligner.ensure_topological().unwrap();

    aligner.refresh_node_indices().unwrap();
    aligner.refresh_node_remaining().unwrap();
}
