#ifndef SGAT_NAVARRO_H
#define SGAT_NAVARRO_H

#include "scoring_schema.h"
#include "sequence_batch.h"
#include "sequence_graph.h"
#include "utils.h"

namespace sgat {

// The class implements a generalized Navarro's algorithm to align query to a
// segment-labeled graph. This class is thread-safe. When we talk about edit
// operations, e.g., substitutions, insertions or deletions, we are using the
// path in the graph as a base. So an insertion means a char in the query would
// be inserted into the path in the graph to match the query. And an insertion
// would consume a char inthe query and would not consume a char of a vertex of
// graph. Size type should be set to one of the signed int types including
// int8_t, int16_t, int32_t, int64_t.
template <class GraphSizeType, class QueryLengthType, class ScoreType>
class NavarroAligner {
 public:
  NavarroAligner() = default;
  ~NavarroAligner() = default;

  QueryLengthType
  AlignUsingLinearGapPenaltyWithNavarroAlgorithmSemiGloballyOnTwoDirections(
      uint32_t sequence_index, const SequenceBatch &sequence_batch,
      const SequenceGraph<GraphSizeType> &sequence_graph,
      const ScoringSchema<ScoreType> &scoring_schema) {
    const QueryLengthType max_cost =
        std::max(std::max(scoring_schema.substitution_penalty,
                          scoring_schema.deletion_penalty),
                 scoring_schema.insertion_penalty);

    const GraphSizeType num_vertices =
        sequence_graph.GetNumVerticesInCompactedGraph();
    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    const char *sequence_bases = sequence_batch.GetSequenceAt(sequence_index);

    std::vector<std::vector<ScoreType>> previous_layer;
    std::vector<std::vector<ScoreType>> current_layer;

    for (GraphSizeType vertex_id = 0; vertex_id < num_vertices; ++vertex_id) {
      GraphSizeType vertex_length =
          sequence_graph.GetVertexLengthInCompactedGraph(vertex_id);
      previous_layer.emplace_back(std::vector<ScoreType>(
          vertex_length, sequence_length * max_cost + 1));
      current_layer.emplace_back(std::vector<ScoreType>(vertex_length, 0));
    }

    GraphSizeType num_propagations = 0;

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      // std::cerr << "layer: " << i << "\n";
      // std::cerr << "previous:\n";
      // PrintLayer(previous_layer);
      // std::cerr << "current:\n";
      // PrintLayer(current_layer);
      std::swap(previous_layer, current_layer);
      ComputeLayerWithNavarroAlgorithm(
          sequence_bases[i], sequence_graph, /*start_vertex=*/num_vertices,
          scoring_schema, previous_layer, num_propagations, current_layer);
    }

    QueryLengthType forward_alignment_cost =
        *std::min_element(current_layer[0].begin(), current_layer[0].end());

    for (GraphSizeType vertex_id = 1; vertex_id < num_vertices; ++vertex_id) {
      QueryLengthType forward_alignment_cost_on_vertex = *std::min_element(
          current_layer[vertex_id].begin(), current_layer[vertex_id].end());
      forward_alignment_cost =
          std::min(forward_alignment_cost, forward_alignment_cost_on_vertex);
    }

    // For reverse complement.
    for (GraphSizeType vertex_id = 0; vertex_id < num_vertices; ++vertex_id) {
      GraphSizeType vertex_length =
          sequence_graph.GetVertexLengthInCompactedGraph(vertex_id);
      previous_layer[vertex_id].assign(vertex_length,
                                       sequence_length * max_cost + 1);
      current_layer[vertex_id].assign(vertex_length, 0);
    }

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      const char complement_base =
          sequence_batch.GetReverseComplementBaseOfSequenceAt(sequence_index,
                                                              i);
      ComputeLayerWithNavarroAlgorithm(
          complement_base, sequence_graph, /*start_vertex=*/num_vertices,
          scoring_schema, previous_layer, num_propagations, current_layer);
    }

    QueryLengthType reverse_complement_alignment_cost =
        *std::min_element(current_layer[0].begin(), current_layer[0].end());

    for (GraphSizeType vertex_id = 1; vertex_id < num_vertices; ++vertex_id) {
      QueryLengthType reverse_complement_alignment_cost_on_vertex =
          *std::min_element(current_layer[vertex_id].begin(),
                            current_layer[vertex_id].end());
      reverse_complement_alignment_cost =
          std::min(reverse_complement_alignment_cost,
                   reverse_complement_alignment_cost_on_vertex);
    }

    const QueryLengthType min_alignment_cost =
        std::min(forward_alignment_cost, reverse_complement_alignment_cost);

    std::cerr << "Sequence length: " << sequence_length
              << ", forward alignment cost:" << forward_alignment_cost
              << ", reverse complement alignment cost:"
              << reverse_complement_alignment_cost
              << ", alignment cost:" << min_alignment_cost
              << ", num propogations: " << num_propagations << std::endl;
    return min_alignment_cost;
  }

  QueryLengthType ForwardExtendUsingLinearGapPenaltyWithNavarroAlgorithm(
      uint32_t sequence_index, const SequenceBatch &sequence_batch,
      const SequenceGraph<GraphSizeType> &sequence_graph,
      GraphSizeType start_vertex,
      const ScoringSchema<ScoreType> &scoring_schema) {
    const QueryLengthType max_cost =
        std::max(std::max(scoring_schema.substitution_penalty,
                          scoring_schema.deletion_penalty),
                 scoring_schema.insertion_penalty);

    const GraphSizeType num_vertices =
        sequence_graph.GetNumVerticesInCompactedGraph();

    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    const char *sequence_bases = sequence_batch.GetSequenceAt(sequence_index);

    std::vector<std::vector<ScoreType>> previous_layer;
    std::vector<std::vector<ScoreType>> current_layer;

    for (GraphSizeType vertex_id = 0; vertex_id < num_vertices; ++vertex_id) {
      const GraphSizeType vertex_length =
          sequence_graph.GetVertexLengthInCompactedGraph(vertex_id);
      previous_layer.emplace_back(std::vector<ScoreType>(
          vertex_length, sequence_length * max_cost + 1));
      current_layer.emplace_back(std::vector<ScoreType>(
          vertex_length, sequence_length * max_cost + 1));
    }

    current_layer[0][0] = 0;

    GraphSizeType num_propagations = 0;

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      ComputeLayerWithNavarroAlgorithm(
          sequence_bases[i], sequence_graph, start_vertex, scoring_schema,
          previous_layer, num_propagations, current_layer);
    }

    QueryLengthType forward_alignment_cost =
        *std::min_element(current_layer[0].begin(), current_layer[0].end());

    for (GraphSizeType vertex_id = 1; vertex_id < num_vertices; ++vertex_id) {
      QueryLengthType forward_alignment_cost_on_vertex = *std::min_element(
          current_layer[vertex_id].begin(), current_layer[vertex_id].end());
      forward_alignment_cost =
          std::min(forward_alignment_cost, forward_alignment_cost_on_vertex);
    }

    std::cerr << "Sequence length: " << sequence_length
              << ", forward alignment cost:" << forward_alignment_cost
              << ", num propogations: " << num_propagations << std::endl;
    return forward_alignment_cost;
  }

 private:
  // Only for debug.
  void PrintLayer(const std::vector<std::vector<ScoreType>> &layer) {
    GraphSizeType vertex_id = 0;
    for (const auto &node_layer : layer) {
      std::cerr << vertex_id << std::endl;
      for (const auto &score : node_layer) {
        std::cerr << score << " ";
      }
      ++vertex_id;
      std::cerr << std::endl;
    }
  }

  void InitializeLayerForGlobalAlignmentExtension(
      char sequence_base, const SequenceGraph<GraphSizeType> &sequence_graph,
      const GraphSizeType num_deletions, const GraphSizeType to,
      const ScoringSchema<ScoreType> &scoring_schema, ScoreType previous_cost,
      GraphSizeType &num_propagations,
      std::vector<std::vector<ScoreType>> &current_layer) {
    num_propagations += 1;
    GraphSizeType vertex_length =
        sequence_graph.GetVertexLengthInCompactedGraph(to);
    for (GraphSizeType vertex_j = 0; vertex_j < vertex_length; ++vertex_j) {
      const QueryLengthType current_num_deletions = num_deletions + vertex_j;
      const bool is_mismatch =
          (sequence_base !=
           sequence_graph.GetVertexLabelInCompactedGraph(to, vertex_j));
      const ScoreType current_cost =
          previous_cost +
          scoring_schema.deletion_penalty * current_num_deletions +
          (is_mismatch ? scoring_schema.substitution_penalty : 0);

      if (current_layer[to][vertex_j] <= current_cost) {
        break;
      }

      current_layer[to][vertex_j] = current_cost;

      if (vertex_j != vertex_length - 1) {
        continue;
      }

      for (const auto &neighbor :
           sequence_graph.GetNeighborsInCompatedGraph(to)) {
        InitializeLayerForGlobalAlignmentExtension(
            sequence_base, sequence_graph, num_deletions + vertex_length,
            neighbor, scoring_schema, previous_cost, num_propagations,
            current_layer);
      }
    }
  }

  void PropagateWithNavarroAlgorithm(
      const SequenceGraph<GraphSizeType> &sequence_graph,
      const GraphSizeType to, const ScoringSchema<ScoreType> &scoring_schema,
      GraphSizeType &num_propagations,
      std::vector<std::vector<ScoreType>> &current_layer) {
    num_propagations += 1;
    GraphSizeType to_vertex_length =
        sequence_graph.GetVertexLengthInCompactedGraph(to);

    // Propagate deletions within the same vertex.
    bool is_last_updated = false;
    for (GraphSizeType vertex_j = 0; vertex_j < to_vertex_length - 1;
         ++vertex_j) {
      const ScoreType new_deletion_cost =
          current_layer[to][vertex_j] + scoring_schema.deletion_penalty;

      if (current_layer[to][vertex_j + 1] > new_deletion_cost) {
        current_layer[to][vertex_j + 1] = new_deletion_cost;
        if (vertex_j + 1 == to_vertex_length - 1) {
          is_last_updated = true;
        }
      } else {
        break;
      }
    }

    if (!is_last_updated) {
      return;
    }

    // Propagate the case of deletions for neighbors.
    for (const auto &neighbor :
         sequence_graph.GetNeighborsInCompatedGraph(to)) {
      const ScoreType new_deletion_cost =
          current_layer[to][to_vertex_length - 1] +
          scoring_schema.deletion_penalty;

      if (current_layer[neighbor][0] > new_deletion_cost) {
        current_layer[neighbor][0] = new_deletion_cost;
        PropagateWithNavarroAlgorithm(sequence_graph, neighbor, scoring_schema,
                                      num_propagations, current_layer);
      }
    }
  }

  void ComputeLayerWithNavarroAlgorithm(
      char sequence_base, const SequenceGraph<GraphSizeType> &sequence_graph,
      GraphSizeType start_vertex,
      const ScoringSchema<ScoreType> &scoring_schema,
      const std::vector<std::vector<ScoreType>> &previous_layer,
      GraphSizeType &num_propagations,
      std::vector<std::vector<ScoreType>> &current_layer) {
    const GraphSizeType num_vertices =
        sequence_graph.GetNumVerticesInCompactedGraph();

    // Initialize current layer
    current_layer[0][0] =
        previous_layer[0][0] + scoring_schema.insertion_penalty;

    for (GraphSizeType vertex_id = 1; vertex_id < num_vertices; ++vertex_id) {
      const GraphSizeType vertex_length =
          sequence_graph.GetVertexLengthInCompactedGraph(vertex_id);

      for (GraphSizeType vertex_j = 0; vertex_j < vertex_length; ++vertex_j) {
        // Handle the case of insertions.
        const ScoreType new_insertion_cost =
            previous_layer[vertex_id][vertex_j] +
            scoring_schema.insertion_penalty;

        if (start_vertex == num_vertices) {
          const bool is_mismatch =
              sequence_base != sequence_graph.GetVertexLabelInCompactedGraph(
                                   vertex_id, vertex_j);

          ScoreType new_start_cost =
              previous_layer[0][0] +
              (is_mismatch ? scoring_schema.substitution_penalty : 0);

          const ScoreType min_cost =
              std::min(new_start_cost, new_insertion_cost);
          current_layer[vertex_id][vertex_j] = min_cost;
        } else {
          current_layer[vertex_id][vertex_j] = new_insertion_cost;
        }
      }
    }

    if (start_vertex != num_vertices) {
      InitializeLayerForGlobalAlignmentExtension(
          sequence_base, sequence_graph,
          /*num_deletions=*/0, start_vertex, scoring_schema,
          previous_layer[0][0], num_propagations, current_layer);
    }

    // std::cerr << "after init:\n";
    // std::cerr << "previous:\n";
    // PrintLayer(previous_layer);
    // std::cerr << "current:\n";
    // PrintLayer(current_layer);

    for (GraphSizeType vertex_id = 1; vertex_id < num_vertices; ++vertex_id) {
      const GraphSizeType vertex_length =
          sequence_graph.GetVertexLengthInCompactedGraph(vertex_id);

      for (GraphSizeType vertex_j = 0; vertex_j < vertex_length - 1;
           ++vertex_j) {
        // Handle the case of substitutions and deletions within the same
        // vertex.
        const bool is_mismatch =
            sequence_base != sequence_graph.GetVertexLabelInCompactedGraph(
                                 vertex_id, vertex_j + 1);

        const ScoreType new_substitution_cost =
            previous_layer[vertex_id][vertex_j] +
            (is_mismatch ? scoring_schema.substitution_penalty : 0);

        const ScoreType new_deletion_cost = current_layer[vertex_id][vertex_j] +
                                            scoring_schema.deletion_penalty;

        const ScoreType min_cost =
            std::min(new_substitution_cost, new_deletion_cost);

        if (current_layer[vertex_id][vertex_j + 1] > min_cost) {
          current_layer[vertex_id][vertex_j + 1] = min_cost;
        }
      }

      // Handle the case of substitutions and deletions for neighbors.
      for (const auto &neighbor :
           sequence_graph.GetNeighborsInCompatedGraph(vertex_id)) {
        const bool is_mismatch =
            sequence_base !=
            sequence_graph.GetVertexLabelInCompactedGraph(neighbor, 0);

        const ScoreType new_substitution_cost =
            previous_layer[vertex_id][vertex_length - 1] +
            (is_mismatch ? scoring_schema.substitution_penalty : 0);

        const ScoreType new_deletion_cost =
            current_layer[vertex_id][vertex_length - 1] +
            scoring_schema.deletion_penalty;

        const ScoreType min_cost =
            std::min(new_substitution_cost, new_deletion_cost);

        if (current_layer[neighbor][0] > min_cost) {
          current_layer[neighbor][0] = min_cost;
        }
      }
    }

    // std::cerr << "before propagate:\n";
    // std::cerr << "previous:\n";
    // PrintLayer(previous_layer);
    // std::cerr << "current:\n";
    // PrintLayer(current_layer);

    for (GraphSizeType vertex_id = 1; vertex_id < num_vertices; ++vertex_id) {
      PropagateWithNavarroAlgorithm(sequence_graph, vertex_id, scoring_schema,
                                    num_propagations, current_layer);
    }
  }
};

}  // namespace sgat
#endif  // SGAT_NAVARRO_H
