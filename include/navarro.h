#ifndef SGAT_NAVARRO_H
#define SGAT_NAVARRO_H

#include "scoring_schema.h"
#include "sequence_batch.h"
#include "sequence_graph.h"
#include "utils.h"

namespace sgat {

template <class GraphSizeType = int32_t, class QueryLengthType = int16_t,
          class ScoreType = int16_t>
class NavarroAligner {
 public:
  NavarroAligner() = default;
  ~NavarroAligner() = default;

  QueryLengthType AlignUsingLinearGapPenaltyWithNavarroAlgorithm(
      uint32_t sequence_index, const SequenceBatch &sequence_batch,
      const SequenceGraph<GraphSizeType> &sequence_graph,
      const ScoringSchema<ScoreType> &scoring_schema) {
    QueryLengthType max_cost =
        std::max(std::max(scoring_schema.substitution_penalty,
                          scoring_schema.deletion_penalty),
                 scoring_schema.insertion_penalty);
    const GraphSizeType num_vertices = sequence_graph.GetNumVertices();
    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    const char *sequence_bases = sequence_batch.GetSequenceAt(sequence_index);
    std::vector<QueryLengthType> previous_layer(num_vertices,
                                                sequence_length * max_cost + 1);
    std::vector<QueryLengthType> current_layer(num_vertices, 0);

    GraphSizeType num_propagations = 0;

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      ComputeLayerWithNavarroAlgorithm(sequence_bases[i], sequence_graph,
                                       scoring_schema, previous_layer,
                                       num_propagations, current_layer);
    }

    const QueryLengthType forward_alignment_cost =
        *std::min_element(current_layer.begin(), current_layer.end());

    // For reverse complement.
    previous_layer.assign(num_vertices, sequence_length * max_cost + 1);
    current_layer.assign(num_vertices, 0);

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      const char complement_base =
          sequence_batch.GetReverseComplementBaseOfSequenceAt(sequence_index,
                                                              i);
      ComputeLayerWithNavarroAlgorithm(complement_base, sequence_graph,
                                       scoring_schema, previous_layer,
                                       num_propagations, current_layer);
    }

    const QueryLengthType reverse_complement_alignment_cost =
        *std::min_element(current_layer.begin(), current_layer.end());

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
    QueryLengthType max_cost =
        std::max(std::max(scoring_schema.substitution_penalty,
                          scoring_schema.deletion_penalty),
                 scoring_schema.insertion_penalty);
    const GraphSizeType num_vertices = sequence_graph.GetNumVertices();
    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    const char *sequence_bases = sequence_batch.GetSequenceAt(sequence_index);
    std::vector<QueryLengthType> previous_layer(num_vertices,
                                                sequence_length * max_cost + 1);
    std::vector<QueryLengthType> current_layer(num_vertices,
                                               sequence_length * max_cost + 1);
    current_layer[0] = scoring_schema.deletion_penalty;
    current_layer[start_vertex] =
        sequence_bases[0] == sequence_graph.GetVertexLabel(start_vertex)
            ? 0
            : scoring_schema.substitution_penalty;

    GraphSizeType num_propagations = 0;

    for (QueryLengthType i = 1; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      ComputeLayerWithNavarroAlgorithm(sequence_bases[i], sequence_graph,
                                       scoring_schema, previous_layer,
                                       num_propagations, current_layer);
    }

    const QueryLengthType forward_alignment_cost =
        *std::min_element(current_layer.begin(), current_layer.end());

    std::cerr << "Sequence length: " << sequence_length
              << ", forward alignment cost:" << forward_alignment_cost
              << ", num propogations: " << num_propagations << std::endl;
    return forward_alignment_cost;
  }

 private:
  void PropagateWithNavarroAlgorithm(
      const SequenceGraph<GraphSizeType> &sequence_graph,
      const GraphSizeType from, const GraphSizeType to,
      const ScoringSchema<ScoreType> &scoring_schema,
      GraphSizeType &num_propagations,
      std::vector<QueryLengthType> &current_layer) {
    num_propagations += 1;
    if (current_layer[to] >
        scoring_schema.insertion_penalty + current_layer[from]) {
      current_layer[to] =
          scoring_schema.insertion_penalty + current_layer[from];
      for (const auto &neighbor : sequence_graph.GetNeighbors(to)) {
        PropagateWithNavarroAlgorithm(sequence_graph, to, neighbor,
                                      scoring_schema, num_propagations,
                                      current_layer);
      }
    }
  }

  void ComputeLayerWithNavarroAlgorithm(
      const char sequence_base,
      const SequenceGraph<GraphSizeType> &sequence_graph,
      const ScoringSchema<ScoreType> &scoring_schema,
      const std::vector<QueryLengthType> &previous_layer,
      GraphSizeType &num_propagations,
      std::vector<QueryLengthType> &current_layer) {
    const GraphSizeType num_vertices = sequence_graph.GetNumVertices();

    // Initialize current layer
    current_layer[0] = previous_layer[0] + scoring_schema.deletion_penalty;
    for (GraphSizeType j = 1; j < num_vertices; ++j) {
      QueryLengthType cost = 0;

      if (sequence_base != sequence_graph.GetVertexLabel(j)) {
        cost = scoring_schema.substitution_penalty;
      }

      current_layer[j] = previous_layer[0] + cost;
    }

    for (GraphSizeType i = 1; i < num_vertices; ++i) {
      if (current_layer[i] >
          previous_layer[i] + scoring_schema.deletion_penalty) {
        current_layer[i] = previous_layer[i] + scoring_schema.deletion_penalty;
      }

      for (const auto &neighbor : sequence_graph.GetNeighbors(i)) {
        QueryLengthType cost = 0;

        if (sequence_base != sequence_graph.GetVertexLabel(neighbor)) {
          cost = scoring_schema.substitution_penalty;
        }

        if (current_layer[neighbor] > previous_layer[i] + cost) {
          current_layer[neighbor] = previous_layer[i] + cost;
        }
      }
    }

    for (GraphSizeType i = 1; i < num_vertices; ++i) {
      for (const auto &neighbor : sequence_graph.GetNeighbors(i)) {
        PropagateWithNavarroAlgorithm(sequence_graph, i, neighbor,
                                      scoring_schema, num_propagations,
                                      current_layer);
      }
    }
  }
};

}  // namespace sgat
#endif  // SGAT_NAVARRO_H
