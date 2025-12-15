//! GraphViz DOT output for an abPOA partial order graph

use crate::encode::{Alphabet, decode_aa_code, decode_dna_code};
use crate::graph::Graph;
use crate::params::SentinelNode;
use crate::{Error, Result};
use std::io::Write;

/// GraphViz rank direction for DOT output
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RankDir {
    LeftToRight,
    TopToBottom,
    RightToLeft,
    BottomToTop,
}

impl RankDir {
    fn as_dot(self) -> &'static str {
        match self {
            RankDir::LeftToRight => "LR",
            RankDir::TopToBottom => "TB",
            RankDir::RightToLeft => "RL",
            RankDir::BottomToTop => "BT",
        }
    }
}

/// How to map abPOA edge weights to GraphViz `penwidth`
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EdgePenWidth {
    /// Mirror upstream `abpoa_dump_pog`: `penwidth = weight + 1`, does not look good for large weights
    Weight,
    /// Scale `weight / max_weight` into `[min, max]`, default
    ProportionalToMax { min: f32, max: f32 },
    /// Scale `weight / sum(outgoing_weights)` into `[min, max]`
    ProportionalToOutgoing { min: f32, max: f32 },
}

impl EdgePenWidth {
    fn compute(self, weight: i32, max_weight: i32, outgoing_sum: i32) -> f32 {
        match self {
            EdgePenWidth::Weight => weight.saturating_add(1).max(1) as f32,
            EdgePenWidth::ProportionalToMax { min, max } => {
                let (min, max) = normalize_range(min, max);
                if max_weight <= 0 {
                    return min;
                }
                let ratio = (weight.max(0) as f32) / (max_weight as f32);
                min + ratio.clamp(0.0, 1.0) * (max - min)
            }
            EdgePenWidth::ProportionalToOutgoing { min, max } => {
                let (min, max) = normalize_range(min, max);
                if outgoing_sum <= 0 {
                    return min;
                }
                let ratio = (weight.max(0) as f32) / (outgoing_sum as f32);
                min + ratio.clamp(0.0, 1.0) * (max - min)
            }
        }
    }
}

/// Which label, if any, to include on each DOT edge
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeLabel {
    None,
    Weight,
    ProportionOfMax { precision: usize },
    ProportionOfOutgoing { precision: usize },
}

/// Which metadata to include in DOT node labels
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeLabel {
    /// Match upstream `abpoa_dump_pog`: base + topological index
    BaseAndTopologicalIndex,
    /// Base + explicit node id (`id:<node_id>`)
    BaseAndNodeId,
    /// Base + explicit topological index (`idx:<topo>`) + explicit node id (`id:<node_id>`)
    BaseTopologicalIndexAndNodeId,
}

/// Rendering settings for [`crate::Aligner::write_pog_to_path`] when using
/// [`crate::PogWriteOptions::Rust`]
#[derive(Debug, Clone, PartialEq)]
pub struct PogDotOptions {
    pub rankdir: RankDir,
    pub node_width: f32,
    pub node_style: &'static str,
    pub node_fixedsize: bool,
    pub node_shape: &'static str,
    pub node_font_size: u32,
    pub node_label: NodeLabel,
    pub edge_font_size: u32,
    pub edge_font_color: &'static str,
    pub edge_label: EdgeLabel,
    pub edge_penwidth: EdgePenWidth,
    pub show_aligned_mismatch: bool,
}

impl Default for PogDotOptions {
    fn default() -> Self {
        Self {
            rankdir: RankDir::LeftToRight,
            node_width: 1.0,
            node_style: "filled",
            node_fixedsize: true,
            node_shape: "circle",
            node_font_size: 24,
            node_label: NodeLabel::BaseAndNodeId,
            edge_font_size: 20,
            edge_font_color: "red",
            edge_label: EdgeLabel::Weight,
            edge_penwidth: EdgePenWidth::ProportionalToMax {
                min: 0.5,
                max: 10.0,
            },
            show_aligned_mismatch: true,
        }
    }
}

/// Write a GraphViz DOT view of a partial order graph (POG)
///
/// This mirrors the upstream `abpoa_dump_pog` DOT emission, but does not invoke `dot`
///
/// Note: the graph must be topologically sorted (e.g. via [`crate::Aligner::ensure_topological`])
pub fn write_pog_dot(
    graph: &Graph<'_>,
    alphabet: Alphabet,
    writer: &mut impl Write,
    options: &PogDotOptions,
) -> Result<()> {
    if graph.is_empty() {
        return Ok(());
    }

    let node_count = graph.node_count();

    let mut max_edge_weight: i32 = 0;
    let mut outgoing_sums = vec![0i32; node_count];
    for topo_idx in 0..node_count {
        let node_id = graph.node_id_at_topological_index(topo_idx)?;
        let node = graph.node(node_id)?;
        let outgoing_sum: i32 = node.out_edges.iter().map(|(_, w)| *w).sum();
        let id_idx = usize::try_from(node_id.0).map_err(|_| {
            Error::InvalidInput("abpoa returned invalid node id in topology buffers".into())
        })?;
        if id_idx >= outgoing_sums.len() {
            return Err(Error::InvalidInput(
                "abpoa returned invalid node id in topology buffers".into(),
            ));
        }
        outgoing_sums[id_idx] = outgoing_sum;
        for (_, weight) in node.out_edges.iter().copied() {
            max_edge_weight = max_edge_weight.max(weight);
        }
    }

    writeln!(writer, "// abpoa graph dot file.\n// {} nodes.", node_count)?;
    writeln!(writer, "digraph ABPOA_graph {{")?;
    writeln!(
        writer,
        "\tgraph [rankdir=\"{}\"];",
        options.rankdir.as_dot()
    )?;
    writeln!(
        writer,
        "\tnode [width={}, style={}, fixedsize={}, shape={}];",
        options.node_width,
        options.node_style,
        if options.node_fixedsize {
            "true"
        } else {
            "false"
        },
        options.node_shape
    )?;

    let source = SentinelNode::Source.as_node_id();
    let sink = SentinelNode::Sink.as_node_id();

    for topo_idx in 0..node_count {
        let node_id = graph.node_id_at_topological_index(topo_idx)?;
        let node = graph.node(node_id)?;
        let (base, color) = if node_id == source {
            ('S', "gray")
        } else if node_id == sink {
            ('E', "gray")
        } else {
            let base = decode_base(alphabet, node.base);
            let color = match alphabet {
                Alphabet::Dna => match node.base {
                    0 => "pink1",     // A
                    1 => "red1",      // C
                    2 => "gold2",     // G
                    3 => "seagreen4", // T
                    _ => "gray",      // N / sentinel / unknown
                },
                Alphabet::AminoAcid => "gray",
            };
            (base, color)
        };

        let label = match options.node_label {
            NodeLabel::BaseAndTopologicalIndex => format!("{base}\\n{topo_idx}"),
            NodeLabel::BaseAndNodeId => format!("{base}\\nid:{}", node_id.0),
            NodeLabel::BaseTopologicalIndexAndNodeId => {
                format!("{base}\\nidx:{topo_idx}\\nid:{}", node_id.0)
            }
        };

        writeln!(
            writer,
            "\tn{} [label=\"{}\", color={}, fontsize={}];",
            node_id.0, label, color, options.node_font_size
        )?;
    }

    let mut x_index: isize = -1;
    for topo_idx in 0..node_count {
        let node_id = graph.node_id_at_topological_index(topo_idx)?;
        let node = graph.node(node_id)?;

        let outgoing_sum = outgoing_sums
            .get(usize::try_from(node_id.0).unwrap_or(usize::MAX))
            .copied()
            .unwrap_or(0);

        for (to, weight) in node.out_edges.iter().copied() {
            let penwidth = options
                .edge_penwidth
                .compute(weight, max_edge_weight, outgoing_sum);
            write!(writer, "\tn{} -> n{} [", node_id.0, to.0)?;

            match options.edge_label {
                EdgeLabel::None => {}
                EdgeLabel::Weight => {
                    write!(writer, "label=\"{}\", ", weight)?;
                }
                EdgeLabel::ProportionOfMax { precision } => {
                    let ratio = if max_edge_weight <= 0 {
                        0.0
                    } else {
                        (weight.max(0) as f32) / (max_edge_weight as f32)
                    };
                    write!(writer, "label=\"{:.*}\", ", precision, ratio)?;
                }
                EdgeLabel::ProportionOfOutgoing { precision } => {
                    let ratio = if outgoing_sum <= 0 {
                        0.0
                    } else {
                        (weight.max(0) as f32) / (outgoing_sum as f32)
                    };
                    write!(writer, "label=\"{:.*}\", ", precision, ratio)?;
                }
            }

            writeln!(
                writer,
                "fontsize={}, fontcolor={}, penwidth={:.3}];",
                options.edge_font_size, options.edge_font_color, penwidth
            )?;
        }

        let aligned = graph.aligned_nodes(node_id)?;
        if !aligned.is_empty() {
            write!(writer, "\t{{rank=same; n{} ", node_id.0)?;
            for other in aligned.iter() {
                write!(writer, "n{} ", other.0)?;
            }
            writeln!(writer, "}};")?;

            if options.show_aligned_mismatch && (topo_idx as isize) > x_index {
                x_index = topo_idx as isize;
                write!(
                    writer,
                    "\t{{ edge [style=dashed, arrowhead=none]; n{} ",
                    node_id.0
                )?;
                for other in aligned.iter() {
                    write!(writer, "-> n{} ", other.0)?;
                    let idx = graph.node_index(*other)?;
                    x_index = x_index.max(idx as isize);
                }
                writeln!(writer, "}}")?;
            }
        }
    }

    writeln!(writer, "}}")?;
    Ok(())
}

fn decode_base(alphabet: Alphabet, code: u8) -> char {
    match alphabet {
        Alphabet::Dna => decode_dna_code(code),
        Alphabet::AminoAcid => decode_aa_code(code),
    }
}

fn normalize_range(min: f32, max: f32) -> (f32, f32) {
    if min.is_finite() && max.is_finite() {
        if min <= max { (min, max) } else { (max, min) }
    } else if min.is_finite() {
        (min, min)
    } else if max.is_finite() {
        (max, max)
    } else {
        (1.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_output_surfaces_consensus_node_ids() {
        let sequences = [b"ACGT".as_ref(), b"ACGT".as_ref()];
        let mut aligner = crate::Aligner::new().unwrap();
        aligner
            .msa_in_place(crate::SequenceBatch::from_sequences(&sequences))
            .unwrap();
        let result = aligner.finalize_msa(crate::OutputMode::CONSENSUS).unwrap();
        assert!(
            !result.clusters.is_empty(),
            "consensus should produce at least one cluster"
        );
        assert!(
            !result.clusters[0].node_ids.is_empty(),
            "consensus should include node ids"
        );

        aligner.ensure_topological().unwrap();
        let graph = aligner.graph().unwrap();

        let mut options = PogDotOptions::default();
        options.node_label = NodeLabel::BaseTopologicalIndexAndNodeId;

        let mut dot = Vec::new();
        write_pog_dot(&graph, crate::Alphabet::Dna, &mut dot, &options).unwrap();
        let dot = String::from_utf8(dot).unwrap();

        assert!(
            dot.contains("\\nidx:") && dot.contains("\\nid:"),
            "node labels should include both topological index and node id"
        );

        for node_id in &result.clusters[0].node_ids {
            assert!(
                dot.contains(&format!("n{} [", node_id.0)),
                "dot output should include consensus node id {}",
                node_id.0
            );
        }

        let first = result.clusters[0].node_ids[0];
        assert!(
            dot.contains(&format!("id:{}", first.0)),
            "dot labels should include the consensus node id"
        );
    }
}
