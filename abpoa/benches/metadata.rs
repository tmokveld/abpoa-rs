use abpoa::{Aligner, OutputMode, Parameters, SequenceBatch};
use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn bench_metadata_overhead(c: &mut Criterion) {
    let seqs = generate_sequences(200, 150);
    let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

    let names: Vec<String> = (0..seq_refs.len())
        .map(|i| format!("read{}", i + 1))
        .collect();
    let name_refs: Vec<&str> = names.iter().map(|s| s.as_str()).collect();

    c.bench_function("msa_no_names", |b| {
        b.iter_batched(
            || {
                let params = Parameters::new().unwrap();
                Aligner::with_params(params).unwrap()
            },
            |mut aligner| {
                aligner
                    .msa(
                        SequenceBatch::from_sequences(&seq_refs),
                        OutputMode::consensus_and_msa(),
                    )
                    .unwrap();
            },
            BatchSize::SmallInput,
        )
    });

    c.bench_function("msa_with_names", |b| {
        b.iter_batched(
            || {
                let params = Parameters::new().unwrap();
                Aligner::with_params(params).unwrap()
            },
            |mut aligner| {
                aligner
                    .msa(
                        SequenceBatch::from_sequences(&seq_refs).with_names(&name_refs),
                        OutputMode::consensus_and_msa(),
                    )
                    .unwrap();
            },
            BatchSize::SmallInput,
        )
    });

    c.bench_function("incremental_add_sequences_with_names", |b| {
        let midpoint = seq_refs.len() / 2;
        let (first, second) = seq_refs.split_at(midpoint);
        let (first_names, second_names) = name_refs.split_at(midpoint);
        b.iter_batched(
            || Aligner::with_params(Parameters::new().unwrap()).unwrap(),
            |mut aligner| {
                aligner
                    .msa_in_place(SequenceBatch::from_sequences(first).with_names(first_names))
                    .unwrap();
                aligner
                    .add_sequences(SequenceBatch::from_sequences(second).with_names(second_names))
                    .unwrap();
                aligner
                    .finalize_msa(OutputMode::consensus_and_msa())
                    .unwrap();
            },
            BatchSize::SmallInput,
        )
    });
}

fn generate_sequences(count: usize, len: usize) -> Vec<Vec<u8>> {
    const BASES: &[u8] = b"ACGT";
    let mut rng = StdRng::seed_from_u64(42u64);
    (0..count)
        .map(|_| {
            (0..len)
                .map(|_| {
                    let idx = rng.random_range(0..BASES.len());
                    BASES[idx]
                })
                .collect()
        })
        .collect()
}

criterion_group!(benches, bench_metadata_overhead);
criterion_main!(benches);
