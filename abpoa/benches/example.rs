use abpoa::{Aligner, OutputMode, Parameters, SequenceBatch};
use criterion::{BatchSize, Criterion, criterion_group, criterion_main};

fn bench_example_dataset(c: &mut Criterion) {
    let seqs: Vec<&[u8]> = vec![
        b"CGTCAATCTATCGAAGCATACGCGGGCAGAGCCGAAGACCTCGGCAATCCA",
        b"CCACGTCAATCTATCGAAGCATACGCGGCAGCCGAACTCGACCTCGGCAATCAC",
        b"CGTCAATCTATCGAAGCATACGCGGCAGAGCCCGGAAGACCTCGGCAATCAC",
        b"CGTCAATGCTAGTCGAAGCAGCTGCGGCAGAGCCGAAGACCTCGGCAATCAC",
        b"CGTCAATCTATCGAAGCATTCTACGCGGCAGAGCCGACCTCGGCAATCAC",
        b"CGTCAATCTAGAAGCATACGCGGCAAGAGCCGAAGACCTCGGCCAATCAC",
        b"CGTCAATCTATCGGTAAAGCATACGCTCTGTAGCCGAAGACCTCGGCAATCAC",
        b"CGTCAATCTATCTTCAAGCATACGCGGCAGAGCCGAAGACCTCGGCAATC",
        b"CGTCAATGGATCGAGTACGCGGCAGAGCCGAAGACCTCGGCAATCAC",
        b"CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT",
    ];
    c.bench_function("example_c_dataset_msa", |b| {
        b.iter_batched(
            || Aligner::with_params(Parameters::new().unwrap()).unwrap(),
            |mut aligner| {
                aligner
                    .msa(
                        SequenceBatch::from_sequences(&seqs),
                        OutputMode::CONSENSUS | OutputMode::MSA,
                    )
                    .unwrap();
            },
            BatchSize::SmallInput,
        );
    });
}

criterion_group!(benches, bench_example_dataset);
criterion_main!(benches);
