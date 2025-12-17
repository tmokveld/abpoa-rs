# abpoa-rs

Rust bindings and wrapper for the [abPOA](https://github.com/yangao07/abPOA) partial order aligner, not related to: https://github.com/broadinstitute/abpoa-rs.

The wrapper lives in the `abpoa` crate, the raw FFI bindings are available as `abpoa::sys` (or via
the `abpoa-ffi-sys` crate directly; `abpoa_ffi_sys` in Rust code).

## Platform support

`abpoa` currently supports Unix-like targets only; upstream abPOA depends on POSIX APIs.

Building requires a C toolchain plus `clang`/`libclang` (for `bindgen`). The build links against
`zlib`, `libm`, and `pthread`.

If you hit build errors about missing `zlib.h` or `-lz`, set `ZLIB_DIR` to a zlib installation
**prefix** (not the `lib/` directory). Expected layout:

- `"$ZLIB_DIR/include/zlib.h"`
- `"$ZLIB_DIR/lib/libz.*"` (or `"$ZLIB_DIR/lib64/libz.*"`)

```sh
export ZLIB_DIR=/path/to/zlib/prefix
```

## Build / test

```sh
cargo build --workspace
cargo test --workspace
```

## Quick start

```rust
use abpoa::{Aligner, Parameters, SequenceBatch};

fn main() -> abpoa::Result<()> {
    let params = Parameters::new()?;

    let mut aligner = Aligner::with_params(params)?;
    let sequences = [b"ACGT".as_ref(), b"ACGG".as_ref(), b"ACGA".as_ref()];
    let result = aligner.msa(SequenceBatch::from_sequences(&sequences))?;
    println!("Consensus: {}", result.clusters[0].consensus);
    Ok(())
}
```

`Parameters` owns all alignment/scoring settings, create it with `Parameters::new()` (defaults)
or `Parameters::configure()` if you want to tweak values before using it. The aligner owns its parameters, configure them up front and pass them into
`Aligner::with_params`.

Minimal configuration example:

```rust
use abpoa::{Aligner, AlignMode, ConsensusAlgorithm, OutputMode, Parameters, Scoring, SequenceBatch};

let seqs = [b"ACGT".as_ref(), b"ACGG".as_ref(), b"ACGA".as_ref()];

let mut params = Parameters::configure()?;
params
    .set_outputs(OutputMode::CONSENSUS)
    .set_align_mode(AlignMode::Global)
    .set_scoring_scheme(Scoring::affine(2, 4, 4, 1))?
    .set_consensus(ConsensusAlgorithm::MostFrequent, 1, 0.0)?;

let mut aligner = Aligner::with_params(params)?;
let result = aligner.msa(SequenceBatch::from_sequences(&seqs))?;
```

## Parameters defaults

`Parameters::new()` / `Parameters::configure()` start from the vendored abPOA defaults with a few
tweaks. The effective defaults are:

- Alphabet: DNA (`Alphabet::Dna`, `m = 5`).
- Outputs: consensus + MSA (`OutputMode::CONSENSUS | OutputMode::MSA`).
- Read ids: enabled (`set_use_read_ids(true)`).
- Alignment mode: global (`AlignMode::Global`).
- Scoring: `Scoring::convex(2, 4, 4, 2, 24, 1)` (match, mismatch, gap_open1, gap_ext1, second_gap_open, second_gap_extend).
- Adaptive band extras: `set_band(10, 0.01)`.
- Z-drop / end bonus: disabled (`None`, stored as `-1` in abPOA).
- Minimizer seeding: disabled (`set_disable_seeding(true)`), with `k=19`, `w=10`, `min_w=500`.
- Progressive POA / guide tree: disabled (`set_progressive_poa(false)`).
- Input sorting: disabled (`set_sort_input_seq(false)`; deprecated/no effect in Rust APIs).
- Ambiguous strand handling: disabled (`set_ambiguous_strand(false)`).
- Gap placement: minimap2-like left-most (`set_gap_on_right(false)`, `set_gap_at_end(false)`).
- Sub-alignment: disabled (`set_sub_alignment(false)`).
- Node-coverage path scoring: disabled (`set_inc_path_score(false)`).
- Quality-weighted scoring: disabled (`set_use_quality(false)`).
- Consensus: `ConsensusAlgorithm::HeaviestBundle`, `max_n_cons = 1`, `min_freq = 0.25`.
- Verbosity: `Verbosity::None`.
- Score matrix file / custom matrix / incremental graph file: unset.

Note: upstream abPOA defaults to MSA output disabled and `use_read_ids = false`. This wrapper
enables both by default so `msa` can return per-read rows, and so MSA/GFA/multi-consensus outputs
can be enabled without rebuilding the graph.

## User tips

Input ordering can affect POA/MSA results, especially with incremental APIs, where progressive POA is not supported
or in low complexity sequences where minimizers may be unreliable. For more consistent
output (and sometimes better alignments), try sorting sequences by decreasing length (longest ->
shortest) before alignment, this is a heuristic used by [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/bar/impl/poaBarAligner.c#L1071-L1072).

## Rust wrapper API

The crate mirrors upstream abPOA, but wraps pointers and lifetimes safely. The main parts are:

- `Aligner`: stateful wrapper around abPOA POA graph. Use it to run alignments, grow/reset the graph,
  and export results.
- `Parameters`: scoring and algorithm configuration (alignment mode, seeding, consensus, RC handling,
  verbosity, etc.), nearly all setters are chainable.
- `SequenceBatch`: a view over input sequences plus optional read names and quality weights.
- `OutputMode`: choose which outputs to compute (`OutputMode::CONSENSUS`, `OutputMode::MSA`, or both) via `Parameters::set_outputs`.
- `MsaResult` / `EncodedMsaResult`: owned alignment output. `msa` holds per-read aligned rows; `clusters`
  holds one or more consensus sequences with coverage/read-id metadata.
- `EncodedMsaView`: zero-copy output view borrowing from an `Aligner`.
- `Graph`: read-only view of the internal POA graph for inspection or subgraph workflows.
- `abpoa::encode` helpers (`encode_dna`, `encode_aa`, `decode_dna`, `decode_aa`) and `Alphabet`: use
  these when working with encoded outputs (`*encoded` APIs) or low-level methods like
  `align_sequence_raw` that take encoded bases.

See further down for examples using the API.

## Thread safety

abPOA is a C library with some global state. The wrapper is conservative about what it allows across
threads and includes a small amount of internal synchronization to avoid known upstream data races.

- `Aligner` and `Parameters` are intentionally `!Send`/`!Sync`: don't share a single instance across
  threads. Create one `Aligner` per worker thread and keep it on that thread.
- Don't run DNA and amino-acid alignments concurrently in the same process. Upstream abPOA uses
  global character mapping tables that are switched during parameter finalization.
- Multi-threading is expected to work well for DNA or amino-acid only workloads when you use one aligner per
  thread (see `abpoa/examples/multi_thread_aligners.rs`).

Ongoing work to improve this.

## One-shot vs incremental APIs

abPOA exposes two different execution paths, and the wrapper mirrors that split:

- **One-shot MSA**: `Aligner::msa` and `Aligner::msa_encoded` call the upstream C driver
  `abpoa_msa`. Each call resets the internal graph and rebuilds it from the provided batch.
  This path is the only one that can use minimizer seeding, guide-tree ordering, and anchored
  window alignment.
- **Incremental / low-level**: APIs such as `Aligner::msa_in_place`, `Aligner::add_sequences`,
  `Aligner::align_sequence_raw`, and `Aligner::add_alignment` align reads directly to the
  current graph in the order you provide them. They preserve graph state between calls and
  are meant for streaming or custom graph edits.

**WARNING:** Because the incremental APIs never go through `abpoa_msa`, they **do not** use minimizer seeding
or guide-tree partitioning even if those parameters are set.

## Incremental API gotchas

The incremental/low-level APIs are intentionally close to the upstream C implementation and assume
that certain inputs are valid. If you violate these invariants, upstream abPOA may crash, call
`exit()`, read/write out of bounds, and/or lead to undefined behavior.

Some gotchas to keep in mind:

- `NodeId` values and `SubgraphRange` bounds (including inputs to `Aligner::subgraph_nodes`) must
  refer to nodes in the *current* graph (never make up ids; discard ids/ranges after
  `Aligner::reset` / `Aligner::restore_graph`).
- Methods that take "encoded" bases/sequences (e.g. `Aligner::add_node`,
  `Aligner::align_sequence_raw`, `Aligner::align_sequence_to_subgraph`,
  `Aligner::add_alignment`, `Aligner::add_subgraph_alignment`) expect abPOA symbol codes for the
  current `Alphabet`; every byte must be in `0..alphabet_size` (use `encode_dna` / `encode_aa` /
  `Alphabet` helpers, not ASCII).
- `RawAlignment` values passed to `add_alignment*` must come from the same `Aligner`/graph state (do
  not reuse alignments across different aligners/graphs or after a reset/restore).
- File-driven and upstream I/O helpers can bypass Rust errors: an unreadable
  `Parameters::set_score_matrix_file` path, malformed `restore_graph` input, or
  `PogWriteOptions::Upstream` output can cause upstream abPOA to call `exit()`. Do not use these
  APIs with untrusted paths/files if process termination is unacceptable.
- More generally, upstream abPOA can call `exit()`/abort on some fatal errors 
  inside the C library. Rust `Result` cannot represent those failures.

## Read IDs (`use_read_ids`) and output constraints

abPOA can optionally track per-read membership in the graph (stored as per-edge bitsets of read ids). This is
required for some outputs and algorithms, but can be disabled to save memory if you only need a single
consensus.

Strict rules enforced by the Rust wrapper:

- Disabling `use_read_ids` is a hard opt-out: the wrapper will return `InvalidInput` if you request any of:
  - `OutputMode::MSA`
  - GFA output (`Aligner::write_gfa*`)
  - `max_consensus > 1` (multi-consensus)
  - `ConsensusAlgorithm::MostFrequent`
- Even if parameters enable read ids, if the current graph was built without them, the wrapper will **reject**
  MSA/GFA/multi-consensus/MostFrequent consensus generation with a clear `InvalidInput` error. You must rebuild
  the graph with read ids enabled from the start.
- Upstream allocates some output-related buffers (notably `node_id_to_msa_rank`) based on the requested
  outputs at graph-build time. The wrapper allocates/resizes these lazily, so it is safe to build a
  consensus-only graph **with read ids enabled** and later request MSA or multi-consensus output.
- Changing `use_read_ids` (or `max_consensus`) via `aligner.params_mut()` does **not** retroactively add read
  ids to an already-built graph.

Memory-saving pattern (single consensus only):

```rust
use abpoa::{Aligner, ConsensusAlgorithm, OutputMode, Parameters, SequenceBatch};

let seqs = [b"ACGT".as_ref(), b"ACGG".as_ref(), b"ACGA".as_ref()];

let mut params = Parameters::configure()?;
params
    .set_outputs(OutputMode::CONSENSUS)
    .set_use_read_ids(false) // allow read-id-less graphs
    .set_consensus(ConsensusAlgorithm::HeaviestBundle, 1, 0.0)?; // keep max_consensus = 1

let mut aligner = Aligner::with_params(params)?;
aligner.msa(SequenceBatch::from_sequences(&seqs))?;
aligner.write_consensus_fasta(&mut std::io::stdout())?;

// These will error because the graph was built without read ids:
assert!(aligner.write_msa_fasta(&mut Vec::new()).is_err());
assert!(aligner.write_gfa(&mut Vec::new()).is_err());
```

## Zero-copy output views (encoded)

If you want to avoid allocating/copying the MSA/consensus output, use the encoded view APIs:

- `Aligner::msa_view_encoded(...) -> EncodedMsaView<'_>`
- `Aligner::finalize_msa_view_encoded(...) -> EncodedMsaView<'_>`

`EncodedMsaView` borrows abPOA’s internal output buffers owned by the aligner; it is invalidated by any
subsequent call that regenerates or clears output (e.g. another `msa*`/`finalize_msa*`, `reset`, or graph
restoration). The borrow checker prevents mutating the aligner while a view exists.

## Seeding, guide trees, and `progressive_poa`

These settings only affect the one-shot path (`msa` / `msa_encoded`) and only when
`set_align_mode` is `Global`:

- `Parameters::set_minimizer_seeding(k, w, min_w)` sets the minimizer parameters and implicitly
  enables seeding (`disable_seeding = false`). This activates abPOA anchored-window alignment
  based on minimizer chains.
- `Parameters::set_disable_seeding(true/false)` gates the anchored-window mode. When `true`,
  no minimizer anchors are used.
- `Parameters::set_progressive_poa(true)` controls whether abPOA builds a minimizer-based guide
  tree to reorder reads before alignment.

The combinations behave as in upstream abPOA:

- `set_disable_seeding(true)` + `set_progressive_poa(false)` (default): plain progressive POA in input
  order.
- `set_disable_seeding(true)` + `set_progressive_poa(true)`: **guide-tree ordering only**. abPOA still
  collects minimizers to build the tree, but aligns full reads (no anchored windows).
- `set_disable_seeding(false)` (via `set_minimizer_seeding(..)` or `set_disable_seeding(false)`):
  anchored-window mode. If `set_progressive_poa(true)` is also set, reads are first reordered by
  the guide tree; if not, they are aligned in input order.

Note: minimizer seeding can behave poorly on very low‑complexity or highly repetitive sequences.
If you see truncated or odd consensus sequences in seeded mode, try disabling seeding or using
larger `k/w` values.

## Examples

Look into the `abpoa/examples/` directory to see examples on how to use the API:

- `example_c.rs`: end-to-end one-shot MSA with minimizer seeding, per-base quality weights, and FASTA/GFA outputs.
- `encoded_views.rs`: zero-copy encoded output views (`EncodedMsaView`, `EncodedMsaRows`, `EncodedClusterView`).
- `oneshot_seeding.rs`: when seeding/guide trees do (and do not) apply, compares one-shot vs incremental paths.
- `multi_thread_aligners.rs`: using one aligner per thread and validating thread safety under TSAN.
- `incremental_msa.rs`: streaming/incremental graph growth (`msa_in_place`, `add_sequences`, `finalize_msa`).
- `manual_graph_build.rs`: low-level alignment loop (`reset`, `align_sequence_raw`, `add_alignment`) and custom scoring/consensus.
- `subgraph_slice_alignment.rs` and `per_read_subalignments.rs`: aligning reads to selected subgraphs/windows.
- `graph_io.rs`: exporting MSA/GFA, restoring graphs, read names, and reverse‑complement aware MSA (`set_ambiguous_strand`).
- `simulate_reads.rs`: generating synthetic reads from two consensuses and recovering both consensus sequences.

The wrapper exposes more than the examples show, feel free to explore.

## Disclaimer

This crate provides Rust bindings to the C library [abPOA](https://github.com/yangao07/abPOA). The bindings in this repository are licensed under the MIT License (see LICENSE). abPOA is licensed under the MIT License; see NOTICE for the upstream copyright and license text.
