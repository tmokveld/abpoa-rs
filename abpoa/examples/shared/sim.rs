use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

const DNA: [u8; 4] = *b"ACGT";

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Complexity {
    High,
    Low,
}

impl Default for Complexity {
    fn default() -> Self {
        Self::High
    }
}

#[derive(Clone, Debug)]
pub struct ConsensusOptions {
    pub len: usize,
    pub complexity: Complexity,
    pub k: usize,
    pub max_kmer_occ: u32,
    pub max_homopolymer_run: usize,
    pub seed: Option<u64>,
    pub max_attempts: usize,
}

#[allow(dead_code)]
impl ConsensusOptions {
    pub fn new(len: usize) -> Self {
        Self {
            len,
            ..Default::default()
        }
    }
}

impl Default for ConsensusOptions {
    fn default() -> Self {
        Self {
            len: 1000,
            complexity: Complexity::High,
            k: 4,
            max_kmer_occ: 2,
            max_homopolymer_run: 3,
            seed: Some(42),
            max_attempts: 128,
        }
    }
}

pub fn make_consensus(mut opts: ConsensusOptions) -> Vec<u8> {
    if opts.len == 0 {
        return Vec::new();
    }
    if opts.k == 0 {
        opts.k = 4;
    }
    if opts.max_attempts == 0 {
        opts.max_attempts = 1;
    }

    let mut seed = opts.seed.unwrap_or_else(rand::random::<u64>);
    for attempt in 0..opts.max_attempts {
        let mut rng = StdRng::seed_from_u64(seed);
        match try_make_consensus(&mut rng, &opts) {
            Ok(seq) => return seq,
            Err(_) => {
                seed = seed
                    .wrapping_add(0x9e3779b97f4a7c15)
                    .wrapping_add(attempt as u64);
            }
        }
    }

    let mut rng = StdRng::seed_from_u64(seed);
    (0..opts.len).map(|_| random_base(&mut rng)).collect()
}

fn try_make_consensus(rng: &mut StdRng, opts: &ConsensusOptions) -> Result<Vec<u8>, ()> {
    let mut seq = Vec::with_capacity(opts.len);

    let use_kmer_counts = matches!(opts.complexity, Complexity::High) && (1..=8).contains(&opts.k);
    let mut kmer_counts = if use_kmer_counts {
        vec![0u32; 1usize << (2 * opts.k)]
    } else {
        Vec::new()
    };

    for _ in 0..opts.len {
        let base = choose_next_base(rng, opts, &seq, &mut kmer_counts)?;
        seq.push(base);
        if use_kmer_counts && seq.len() >= opts.k {
            let idx = kmer_index(&seq[seq.len() - opts.k..])?;
            kmer_counts[idx] = kmer_counts[idx].saturating_add(1);
        }
    }

    Ok(seq)
}

fn choose_next_base(
    rng: &mut StdRng,
    opts: &ConsensusOptions,
    prefix: &[u8],
    kmer_counts: &mut [u32],
) -> Result<u8, ()> {
    let mut candidates = DNA;
    for i in (1..candidates.len()).rev() {
        let j = rng.random_range(0..=i);
        candidates.swap(i, j);
    }

    for &b in &candidates {
        if opts.max_homopolymer_run > 0
            && would_violate_homopolymer(prefix, b, opts.max_homopolymer_run)
        {
            continue;
        }
        if matches!(opts.complexity, Complexity::High) {
            if would_violate_kmer(prefix, b, opts, kmer_counts) {
                continue;
            }
        }
        return Ok(b);
    }

    Err(())
}

fn would_violate_homopolymer(prefix: &[u8], next: u8, max_run: usize) -> bool {
    if max_run == 0 {
        return false;
    }
    let mut run = 1usize;
    for &b in prefix.iter().rev() {
        if b == next {
            run += 1;
            if run > max_run {
                return true;
            }
        } else {
            break;
        }
    }
    false
}

fn would_violate_kmer(
    prefix: &[u8],
    next: u8,
    opts: &ConsensusOptions,
    kmer_counts: &[u32],
) -> bool {
    if !(1..=8).contains(&opts.k) {
        return false;
    }
    if prefix.len() + 1 < opts.k {
        return false;
    }

    let start = prefix.len() + 1 - opts.k;
    let mut kmer = Vec::with_capacity(opts.k);
    kmer.extend_from_slice(&prefix[start..]);
    kmer.push(next);

    let Ok(idx) = kmer_index(&kmer) else {
        return false;
    };
    kmer_counts.get(idx).copied().unwrap_or(0) >= opts.max_kmer_occ
}

fn kmer_index(kmer: &[u8]) -> Result<usize, ()> {
    let mut idx: usize = 0;
    for &b in kmer {
        idx <<= 2;
        idx |= match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return Err(()),
        };
    }
    Ok(idx)
}

fn random_base(rng: &mut StdRng) -> u8 {
    DNA[rng.random_range(0..DNA.len())]
}

#[derive(Clone, Copy, Debug)]
pub struct ErrorRates {
    pub sub: f64,
    pub ins: f64,
    pub del: f64,
}

impl ErrorRates {
    pub fn new(sub: f64, ins: f64, del: f64) -> Self {
        Self { sub, ins, del }
    }
}

impl Default for ErrorRates {
    fn default() -> Self {
        Self {
            sub: 0.01,
            ins: 0.0025,
            del: 0.0025,
        }
    }
}

pub struct ReadSimulator {
    consensus: Vec<u8>,
    remaining: usize,
    rates: ErrorRates,
    rng: StdRng,
}

pub fn simulate_reads(
    consensus: &[u8],
    n: usize,
    rates: ErrorRates,
    seed: Option<u64>,
) -> ReadSimulator {
    let seed = seed.unwrap_or_else(rand::random::<u64>);
    ReadSimulator {
        consensus: consensus.to_vec(),
        remaining: n,
        rates,
        rng: StdRng::seed_from_u64(seed),
    }
}

impl Iterator for ReadSimulator {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }
        self.remaining -= 1;

        let mut out = Vec::with_capacity(self.consensus.len());
        for &ref_base in &self.consensus {
            while self.rng.random::<f64>() < self.rates.ins {
                out.push(random_base(&mut self.rng));
            }

            if self.rng.random::<f64>() < self.rates.del {
                continue;
            }

            let mut base = ref_base;
            if self.rng.random::<f64>() < self.rates.sub {
                base = mutate_base(&mut self.rng, ref_base);
            }
            out.push(base);
        }

        while self.rng.random::<f64>() < self.rates.ins {
            out.push(random_base(&mut self.rng));
        }

        Some(out)
    }
}

fn mutate_base(rng: &mut StdRng, base: u8) -> u8 {
    let mut choices = [b'A', b'C', b'G', b'T'];
    let mut n = 0usize;
    for &b in &DNA {
        if b != base {
            choices[n] = b;
            n += 1;
        }
    }
    choices[rng.random_range(0..3)]
}
