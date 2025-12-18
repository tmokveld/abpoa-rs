//! Read-only view of the underlying partial order graph.

use crate::params::NodeId;
use crate::{Error, Result, sys};
use std::{fmt, marker::PhantomData, ptr::NonNull, rc::Rc, slice};

/// View of the metadata and connectivity of a single node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NodeRef {
    pub id: NodeId,
    pub base: u8,
    pub in_edges: Vec<(NodeId, i32)>,
    pub out_edges: Vec<(NodeId, i32)>,
    pub n_read: i32,
    pub n_span_read: i32,
}

/// Borrowed view of the metadata and connectivity of a single node.
///
/// Unlike [`NodeRef`], this view borrows edge arrays from the underlying graph without
/// allocating.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NodeView<'a> {
    pub id: NodeId,
    pub base: u8,
    pub in_ids: &'a [i32],
    pub in_w: &'a [i32],
    pub out_ids: &'a [i32],
    pub out_w: &'a [i32],
    pub n_read: i32,
    pub n_span_read: i32,
}

/// Borrowed view into `abpoa_graph_t` tied to an `Aligner` borrow.
pub struct Graph<'a> {
    graph: NonNull<sys::abpoa_graph_t>,
    seqs: NonNull<sys::abpoa_seq_t>,
    _owner: PhantomData<&'a sys::abpoa_t>,
    // Not Send/Sync: abPOA relies on global mutable tables without synchronization.
    _not_send_sync: PhantomData<Rc<()>>,
}

impl<'a> Graph<'a> {
    pub(crate) fn new(graph: *mut sys::abpoa_graph_t, seqs: *mut sys::abpoa_seq_t) -> Result<Self> {
        let graph =
            NonNull::new(graph).ok_or(Error::NullPointer("abpoa graph pointer was null"))?;
        let seqs =
            NonNull::new(seqs).ok_or(Error::NullPointer("abpoa sequence pointer was null"))?;
        Ok(Self {
            graph,
            seqs,
            _owner: PhantomData,
            _not_send_sync: PhantomData,
        })
    }

    /// Number of nodes present in the graph (+source/sink sentinels).
    pub fn node_count(&self) -> usize {
        // Safety: the graph pointer is valid for the lifetime of this view.
        let n = unsafe { (*self.graph.as_ptr()).node_n };
        n.max(0) as usize
    }

    /// Whether the graph has only sentinel nodes (+no aligned bases).
    pub fn is_empty(&self) -> bool {
        self.node_count() <= 2
    }

    /// Number of sequences recorded in the graph.
    pub fn sequence_count(&self) -> usize {
        // Safety: `n_seq` is owned by the aligner and valid for the lifetime of `self`.
        let n = unsafe { (*self.seqs.as_ptr()).n_seq };
        n.max(0) as usize
    }

    /// Whether the graph has been finalized with a consensus sequence.
    pub fn has_consensus(&self) -> bool {
        unsafe { (*self.graph.as_ptr()).is_called_cons() > 0 }
    }

    /// Iterate over all nodes in the graph (including source and sink).
    pub fn nodes(&self) -> GraphNodes<'_> {
        GraphNodes {
            graph: self,
            next: 0,
            total: self.node_count(),
        }
    }

    /// Node ID at a given topological index.
    ///
    /// Requires the graph to be topologically sorted (e.g. via
    /// [`crate::Aligner::ensure_topological`]).
    pub fn node_id_at_topological_index(&self, index: usize) -> Result<NodeId> {
        let graph = unsafe { self.graph.as_ref() };
        if index >= graph.node_n.max(0) as usize {
            return Err(Error::InvalidInput(
                "topological index out of bounds for current graph".into(),
            ));
        }
        let ptr = graph.index_to_node_id;
        if ptr.is_null() {
            return Err(Error::NullPointer(
                "index_to_node_id was null; call ensure_topological first",
            ));
        }
        // Safety: `index_to_node_id` is allocated by abPOA and sized to `node_n` after
        // topological sorting; `index` is validated against `node_n` above.
        Ok(NodeId(unsafe { *ptr.add(index) }))
    }

    /// Snapshot metadata for a node by its id.
    ///
    /// This allocates edge vectors; for a borrowed, zero-allocation view use.
    /// [`Graph::node_view`].
    pub fn node(&self, id: NodeId) -> Result<NodeRef> {
        let node = self.node_by_id(id)?;
        Ok(self.build_node_ref(node))
    }

    /// Borrow metadata and edge arrays for a node by its id.
    pub fn node_view(&self, id: NodeId) -> Result<NodeView<'_>> {
        let node = self.node_by_id(id)?;
        Ok(self.build_node_view(node))
    }

    /// Sum of outgoing edge weights for the given node.
    pub fn node_weight(&self, id: NodeId) -> Result<i32> {
        let node = self.node_by_id(id)?;
        let weights = self.edge_slice(node.out_edge_weight, node.out_edge_n);
        Ok(weights.iter().copied().sum())
    }

    /// Node ids aligned to the given node (mismatches at the same rank).
    pub fn aligned_nodes(&self, id: NodeId) -> Result<Vec<NodeId>> {
        let node = self.node_by_id(id)?;
        let aligned = self.edge_slice(node.aligned_node_id, node.aligned_node_n);
        Ok(aligned.iter().copied().map(NodeId).collect())
    }

    /// Sequence metadata loaded into the graph.
    pub fn sequences(&self) -> Sequences<'_> {
        // Safety: the sequence container is allocated alongside the graph and lives for `'self`.
        let seqs = unsafe { self.seqs.as_ref() };
        Sequences::new(seqs)
    }

    /// Topological index for a node ID, if the graph has been ordered.
    pub fn node_index(&self, id: NodeId) -> Result<usize> {
        self.read_mapping(
            self.raw_mapping(|g| g.node_id_to_index),
            id,
            "node_id_to_index was null; finalize the graph before reading indices",
        )
    }

    /// Maximum reachable position to the left of the node in topological order.
    pub fn node_max_pos_left(&self, id: NodeId) -> Result<usize> {
        self.read_mapping(
            self.raw_mapping(|g| g.node_id_to_max_pos_left),
            id,
            "node_id_to_max_pos_left was null; finalize the graph before reading positions",
        )
    }

    /// Maximum reachable position to the right of the node in topological order.
    pub fn node_max_pos_right(&self, id: NodeId) -> Result<usize> {
        self.read_mapping(
            self.raw_mapping(|g| g.node_id_to_max_pos_right),
            id,
            "node_id_to_max_pos_right was null; finalize the graph before reading positions",
        )
    }

    /// Node metadata by topological index.
    fn node_by_index(&self, index: usize) -> Option<NodeRef> {
        let graph = unsafe { self.graph.as_ref() };
        if graph.node.is_null() || index >= graph.node_n.max(0) as usize {
            return None;
        }
        let node = unsafe { graph.node.add(index).as_ref()? };
        Some(self.build_node_ref(node))
    }

    /// Node metadata by ID.
    fn node_by_id(&self, id: NodeId) -> Result<&sys::abpoa_node_t> {
        let index = self.validate_node_id(id)?;
        let graph = unsafe { self.graph.as_ref() };
        let node_ptr = graph.node;
        if node_ptr.is_null() {
            return Err(Error::NullPointer("abpoa graph node array was null"));
        }
        // Safety: the node array is owned by the aligner; we validated the index above.
        let node = unsafe { node_ptr.add(index).as_ref() }
            .ok_or(Error::NullPointer("abpoa graph node pointer was null"))?;
        Ok(node)
    }

    /// Build a `NodeRef` from a raw `abpoa_node_t`.
    fn build_node_ref(&self, node: &sys::abpoa_node_t) -> NodeRef {
        NodeRef {
            id: NodeId(node.node_id),
            base: node.base,
            in_edges: self.collect_edges(node.in_id, node.in_edge_weight, node.in_edge_n),
            out_edges: self.collect_edges(node.out_id, node.out_edge_weight, node.out_edge_n),
            n_read: node.n_read.max(0),
            n_span_read: node.n_span_read.max(0),
        }
    }

    /// Build a `NodeView` from a raw `abpoa_node_t`.
    fn build_node_view(&self, node: &sys::abpoa_node_t) -> NodeView<'_> {
        let (in_ids, in_w) = self.edge_slices(node.in_id, node.in_edge_weight, node.in_edge_n);
        let (out_ids, out_w) = self.edge_slices(node.out_id, node.out_edge_weight, node.out_edge_n);
        NodeView {
            id: NodeId(node.node_id),
            base: node.base,
            in_ids,
            in_w,
            out_ids,
            out_w,
            n_read: node.n_read.max(0),
            n_span_read: node.n_span_read.max(0),
        }
    }

    fn validate_node_id(&self, id: NodeId) -> Result<usize> {
        if id.0 < 0 {
            return Err(Error::InvalidInput("node id cannot be negative".into()));
        }
        let count = self.node_count();
        if id.0 as usize >= count {
            return Err(Error::InvalidInput(
                "node id out of bounds for current graph".into(),
            ));
        }
        Ok(id.0 as usize)
    }

    fn raw_mapping(&self, selector: fn(&sys::abpoa_graph_t) -> *mut i32) -> Option<*mut i32> {
        let graph = unsafe { self.graph.as_ref() };
        let ptr = selector(graph);
        if ptr.is_null() { None } else { Some(ptr) }
    }

    fn read_mapping(
        &self,
        mapping: Option<*mut i32>,
        id: NodeId,
        null_msg: &'static str,
    ) -> Result<usize> {
        let index = self.validate_node_id(id)?;
        let ptr = mapping.ok_or(Error::NullPointer(null_msg))?;
        // Safety: mapping arrays are allocated by abPOA and sized to node count; index was
        // validated against node_n above.
        let value = unsafe { *ptr.add(index) };
        if value < 0 {
            return Err(Error::InvalidInput(
                "mapping contained negative value".into(),
            ));
        }
        Ok(value as usize)
    }

    fn edge_slice(&self, ptr: *mut i32, len: i32) -> &[i32] {
        if len <= 0 || ptr.is_null() {
            &[]
        } else {
            // Safety: abPOA allocates edge arrays to at least `len` when len > 0
            unsafe { slice::from_raw_parts(ptr, len as usize) }
        }
    }

    fn edge_slices(&self, ids: *mut i32, weights: *mut i32, count: i32) -> (&[i32], &[i32]) {
        let ids = self.edge_slice(ids, count);
        let weights = self.edge_slice(weights, count);
        let len = ids.len().min(weights.len());
        (&ids[..len], &weights[..len])
    }

    fn collect_edges(&self, ids: *mut i32, weights: *mut i32, count: i32) -> Vec<(NodeId, i32)> {
        let ids = self.edge_slice(ids, count);
        let weights = self.edge_slice(weights, count);
        ids.iter()
            .zip(weights.iter().copied())
            .map(|(&id, weight)| (NodeId(id), weight))
            .collect()
    }
}

/// Collection of sequences attached to the graph.
pub struct Sequences<'a> {
    bases: &'a [sys::abpoa_str_t],
    names: &'a [sys::abpoa_str_t],
    comments: &'a [sys::abpoa_str_t],
    qualities: &'a [sys::abpoa_str_t],
    is_reverse_complement: &'a [u8],
    count: usize,
}

impl<'a> Sequences<'a> {
    fn new(raw: &'a sys::abpoa_seq_t) -> Self {
        let count = raw.n_seq.max(0) as usize;
        Self {
            bases: Self::str_slice(raw.seq, count),
            names: Self::str_slice(raw.name, count),
            comments: Self::str_slice(raw.comment, count),
            qualities: Self::str_slice(raw.qual, count),
            is_reverse_complement: Self::flag_slice(raw.is_rc, count),
            count,
        }
    }

    /// Number of sequences recorded in the graph.
    pub fn len(&self) -> usize {
        self.count
    }

    /// Whether any sequences were provided to the aligner.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Borrow an entry by index.
    pub fn get(&self, index: usize) -> Option<Sequence<'a>> {
        if index < self.len() {
            Some(self.sequence_at(index))
        } else {
            None
        }
    }

    /// Iterate over all sequences and their metadata.
    pub fn iter(&self) -> SequenceIter<'_> {
        SequenceIter {
            sequences: self,
            next: 0,
        }
    }

    /// Iterate over the raw sequence strings.
    pub fn sequences(&self) -> impl Iterator<Item = SequenceStr<'a>> + 'a {
        let bases = self.bases;
        let count = self.count;
        (0..count).map(move |index| SequenceStr::from_raw(bases.get(index)))
    }

    /// Iterate over sequence names (may be empty strings).
    pub fn names(&self) -> impl Iterator<Item = SequenceStr<'a>> + 'a {
        let names = self.names;
        let count = self.count;
        (0..count).map(move |index| SequenceStr::from_raw(names.get(index)))
    }

    /// Iterate over optional sequence comments.
    pub fn comments(&self) -> impl Iterator<Item = SequenceStr<'a>> + 'a {
        let comments = self.comments;
        let count = self.count;
        (0..count).map(move |index| SequenceStr::from_raw(comments.get(index)))
    }

    /// Iterate over qualities if present.
    pub fn qualities(&self) -> impl Iterator<Item = SequenceStr<'a>> + 'a {
        let qualities = self.qualities;
        let count = self.count;
        (0..count).map(move |index| SequenceStr::from_raw(qualities.get(index)))
    }

    /// Iterate over reverse-complement flags in input order.
    pub fn reverse_complements(&self) -> impl Iterator<Item = bool> + 'a {
        self.is_reverse_complement.iter().map(|flag| *flag != 0)
    }

    fn sequence_at(&self, index: usize) -> Sequence<'a> {
        let sequence = SequenceStr::from_raw(self.bases.get(index));
        let name = SequenceStr::from_raw(self.names.get(index));
        let comment = SequenceStr::from_raw(self.comments.get(index));
        let quality = SequenceStr::from_raw(self.qualities.get(index));

        Sequence {
            sequence,
            name,
            comment,
            quality,
            is_reverse_complement: self
                .is_reverse_complement
                .get(index)
                .copied()
                .unwrap_or_default()
                != 0,
        }
    }

    fn str_slice(ptr: *mut sys::abpoa_str_t, len: usize) -> &'a [sys::abpoa_str_t] {
        if len == 0 || ptr.is_null() {
            &[]
        } else {
            // Safety: abPOA allocates `len` entries when `n_seq` is set to at least `len`.
            unsafe { slice::from_raw_parts(ptr, len) }
        }
    }

    fn flag_slice(ptr: *mut u8, len: usize) -> &'a [u8] {
        if len == 0 || ptr.is_null() {
            &[]
        } else {
            // Safety: `is_rc` is sized to `n_seq` entries alongside other sequence buffers.
            unsafe { slice::from_raw_parts(ptr, len) }
        }
    }
}

/// Borrowed view of a single sequence's metadata.
#[derive(Clone, Copy, Debug)]
pub struct Sequence<'a> {
    pub sequence: SequenceStr<'a>,
    pub name: SequenceStr<'a>,
    pub comment: SequenceStr<'a>,
    pub quality: SequenceStr<'a>,
    pub is_reverse_complement: bool,
}

/// Non-owning string wrapper backed by `abpoa_str_t`.
#[derive(Clone, Copy)]
pub struct SequenceStr<'a> {
    raw: Option<NonNull<sys::abpoa_str_t>>,
    _owner: PhantomData<&'a sys::abpoa_seq_t>,
}

impl<'a> SequenceStr<'a> {
    fn from_raw(raw: Option<&'a sys::abpoa_str_t>) -> Self {
        Self {
            raw: raw.map(NonNull::from),
            _owner: PhantomData,
        }
    }

    /// Length in bytes of the underlying buffer.
    pub fn len(&self) -> usize {
        let Some(raw) = self.raw else { return 0 };
        // Safety: pointer was derived from an abPOA-managed slice.
        let raw = unsafe { raw.as_ref() };
        if raw.s.is_null() {
            0
        } else {
            raw.l.max(0) as usize
        }
    }

    /// Whether the buffer is empty or missing.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Raw bytes of the string (may not be valid UTF-8).
    pub fn as_bytes(&self) -> &'a [u8] {
        let Some(raw) = self.raw else {
            return &[];
        };
        // Safety: pointer was derived from an abPOA-managed slice.
        let raw = unsafe { raw.as_ref() };
        let len = raw.l.max(0) as usize;
        if len == 0 || raw.s.is_null() {
            return &[];
        }
        // Safety: abPOA allocates `s` to at least `len` bytes when `l > 0`.
        unsafe { slice::from_raw_parts(raw.s as *const u8, len) }
    }

    /// Interpret the buffer as UTF-8 when possible.
    pub fn as_str(&self) -> Option<&'a str> {
        std::str::from_utf8(self.as_bytes()).ok()
    }
}

impl<'a> AsRef<[u8]> for SequenceStr<'a> {
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl<'a> fmt::Debug for SequenceStr<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(text) = self.as_str() {
            write!(f, "{text:?}")
        } else {
            write!(f, "{:?}", self.as_bytes())
        }
    }
}

/// Iterator over sequences produced by [`Graph::sequences`].
pub struct SequenceIter<'a> {
    sequences: &'a Sequences<'a>,
    next: usize,
}

impl<'a> Iterator for SequenceIter<'a> {
    type Item = Sequence<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let idx = self.next;
        if idx < self.sequences.len() {
            self.next += 1;
            self.sequences.get(idx)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.sequences.len().saturating_sub(self.next);
        (remaining, Some(remaining))
    }
}

impl<'a> ExactSizeIterator for SequenceIter<'a> {
    fn len(&self) -> usize {
        self.sequences.len().saturating_sub(self.next)
    }
}

/// Iterator over graph nodes produced by [`Graph::nodes`].
pub struct GraphNodes<'a> {
    graph: &'a Graph<'a>,
    next: usize,
    total: usize,
}

impl<'a> Iterator for GraphNodes<'a> {
    type Item = NodeRef;

    fn next(&mut self) -> Option<Self::Item> {
        while self.next < self.total {
            let idx = self.next;
            self.next += 1;
            if let Some(node) = self.graph.node_by_index(idx) {
                return Some(node);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aligner::{Aligner, SequenceBatch};
    use crate::params::{Parameters, SentinelNode};

    fn with_graph<F>(seqs: &[&[u8]], f: F)
    where
        F: for<'g> FnOnce(&'g Graph<'g>, Vec<NodeRef>),
    {
        let params = Parameters::new().unwrap();
        let mut aligner = Aligner::with_params(params).unwrap();
        aligner
            .msa_in_place(SequenceBatch::from_sequences(seqs).unwrap())
            .unwrap();
        aligner.finalize_msa().unwrap();
        let graph = aligner.graph().unwrap();
        let nodes: Vec<_> = graph.nodes().collect();
        f(&graph, nodes);
    }

    #[test]
    fn iterates_nodes_and_edges() {
        with_graph(&[b"ACGT", b"ACGA"], |graph, nodes| {
            assert_eq!(nodes.len(), graph.node_count());

            let src_id = SentinelNode::Source.as_node_id();
            let sink_id = SentinelNode::Sink.as_node_id();
            let src = nodes.iter().find(|node| node.id == src_id).unwrap();
            assert!(
                !src.out_edges.is_empty(),
                "source should have outgoing edges once sequences are added"
            );

            let sink = nodes.iter().find(|node| node.id == sink_id).unwrap();
            assert!(
                sink.out_edges.is_empty(),
                "sink should have no outgoing edges"
            );
            assert!(
                !sink.in_edges.is_empty(),
                "sink should have incoming edges from terminal bases"
            );
        });
    }

    #[test]
    fn node_weight_matches_edge_sum() {
        with_graph(&[b"ACGT", b"ACGA"], |graph, nodes| {
            let src_id = SentinelNode::Source.as_node_id();
            let src = nodes.iter().find(|node| node.id == src_id).unwrap();
            let expected: i32 = src.out_edges.iter().map(|(_, w)| *w).sum();
            let weight = graph.node_weight(src_id).unwrap();
            assert_eq!(weight, expected);
        });
    }

    #[test]
    fn aligned_nodes_include_tail_variants() {
        with_graph(&[b"ACGT", b"ACGA"], |graph, nodes| {
            let sink_id = SentinelNode::Sink.as_node_id();
            let src_id = SentinelNode::Source.as_node_id();
            let terminal: Vec<_> = nodes
                .iter()
                .filter(|node| node.id != sink_id && node.id != src_id)
                .filter(|node| node.out_edges.iter().any(|(id, _)| *id == sink_id))
                .collect();
            assert!(
                terminal.len() >= 2,
                "variant tail nodes should exist for mismatching bases"
            );
            let aligned = graph.aligned_nodes(terminal[0].id).unwrap();
            assert!(
                aligned.contains(&terminal[1].id),
                "terminal nodes should be aligned to each other"
            );
        });
    }

    #[test]
    fn node_indices_and_positions_accessible() {
        with_graph(&[b"ACGT", b"ACGA"], |graph, nodes| {
            let src_id = SentinelNode::Source.as_node_id();
            let sink_id = SentinelNode::Sink.as_node_id();
            let src_index = graph.node_index(src_id).unwrap();
            let sink_index = graph.node_index(sink_id).unwrap();
            assert!(src_index < sink_index);

            let first_base = nodes
                .iter()
                .find(|node| node.id != src_id && node.id != sink_id)
                .unwrap();
            let limit = graph.node_count();
            let left = graph.node_max_pos_left(first_base.id).unwrap();
            let right = graph.node_max_pos_right(first_base.id).unwrap();
            assert!(left <= limit);
            assert!(right <= limit);
            assert!(src_index <= limit);
            assert!(sink_index <= limit);
        });
    }

    #[test]
    fn sequences_view_surfaces_names_and_flags() {
        let mut params = Parameters::configure().unwrap();
        params.set_ambiguous_strand(true);
        let mut aligner = Aligner::with_params(params).unwrap();
        let seqs = [b"ACGTTGC".as_ref(), b"GCAACGT".as_ref()];
        let names = ["read1", "read2"];

        aligner
            .msa(SequenceBatch::from_sequences(&seqs).unwrap().with_names(&names).unwrap())
            .unwrap();

        let graph = aligner.graph().unwrap();
        let sequences = graph.sequences();
        assert_eq!(sequences.len(), seqs.len());

        let entries: Vec<_> = sequences.iter().collect();
        assert_eq!(entries[0].name.as_str(), Some(names[0]));
        assert_eq!(entries[1].name.as_str(), Some(names[1]));
        assert!(
            entries.iter().any(|seq| seq.is_reverse_complement),
            "reverse complement flags should be preserved when ambiguous strand is enabled"
        );
    }

    #[test]
    fn sequence_count_matches_sequences_len() {
        with_graph(&[b"ACGT", b"ACGA"], |graph, _nodes| {
            let sequences = graph.sequences();
            assert_eq!(graph.sequence_count(), sequences.len());
            assert_eq!(graph.sequence_count(), sequences.iter().count());
        });
    }

    #[test]
    fn node_by_index_round_trips_node_ids() {
        with_graph(&[b"ACGT", b"ACGA"], |graph, _nodes| {
            let count = graph.node_count().max(0);
            for idx in 0..count {
                let by_index = graph.node_by_index(idx).unwrap();
                let by_id = graph.node(by_index.id).unwrap();
                assert_eq!(by_id.id, by_index.id);
            }
        });
    }

    #[test]
    fn node_view_matches_node_snapshot() {
        with_graph(&[b"ACGT", b"ACGA"], |graph, nodes| {
            for node in nodes {
                let view = graph.node_view(node.id).unwrap();
                assert_eq!(view.id, node.id);
                assert_eq!(view.base, node.base);
                assert_eq!(view.n_read, node.n_read);
                assert_eq!(view.n_span_read, node.n_span_read);
                assert_eq!(view.in_ids.len(), view.in_w.len());
                assert_eq!(view.out_ids.len(), view.out_w.len());

                let view_in: Vec<_> = view
                    .in_ids
                    .iter()
                    .copied()
                    .zip(view.in_w.iter().copied())
                    .collect();
                let expected_in: Vec<_> = node.in_edges.iter().map(|(id, w)| (id.0, *w)).collect();
                assert_eq!(view_in, expected_in);

                let view_out: Vec<_> = view
                    .out_ids
                    .iter()
                    .copied()
                    .zip(view.out_w.iter().copied())
                    .collect();
                let expected_out: Vec<_> =
                    node.out_edges.iter().map(|(id, w)| (id.0, *w)).collect();
                assert_eq!(view_out, expected_out);
            }

            assert!(matches!(
                graph.node_view(NodeId(-1)),
                Err(Error::InvalidInput(msg)) if msg == "node id cannot be negative"
            ));
            let oob_id = NodeId(i32::try_from(graph.node_count()).unwrap());
            assert!(matches!(
                graph.node_view(oob_id),
                Err(Error::InvalidInput(msg)) if msg == "node id out of bounds for current graph"
            ));
        });
    }
}
