#ifndef SGAT_SGA_H
#define SGAT_SGA_H

#include "scoring_schema.h"
#include "sequence_batch.h"
#include "sequence_graph.h"
#include "utils.h"

namespace sgat {

template <class GraphSizeType = int32_t, class QueryLengthType = int16_t,
          class ScoreType = int16_t>
class SgaAligner {
 public:
  SgaAligner() = default;
  ~SgaAligner() = default;

  ScoreType AlignUsingLinearGapPenalty(
      uint32_t sequence_index, const SequenceBatch &sequence_batch,
      const SequenceGraph<GraphSizeType> &sequence_graph,
      const ScoringSchema<ScoreType> &scoring_schema) {
    const GraphSizeType num_vertices = sequence_graph.GetNumVertices();

    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    const char *sequence_bases = sequence_batch.GetSequenceAt(sequence_index);

    const ScoreType max_cost =
        std::max(std::max(scoring_schema.substitution_penalty,
                          scoring_schema.deletion_penalty),
                 scoring_schema.insertion_penalty);

    std::vector<ScoreType> previous_layer(num_vertices);
    std::vector<GraphSizeType> previous_order(num_vertices);
    std::vector<ScoreType> initialized_layer(num_vertices,
                                             sequence_length * max_cost + 1);
    std::vector<GraphSizeType> initialized_order(num_vertices);
    std::vector<ScoreType> current_layer(num_vertices, 0);
    std::vector<GraphSizeType> current_order;
    current_order.reserve(num_vertices);

    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      current_order.push_back(i);
    }

    order_look_up_table_.assign(3 * num_vertices, 0);
    parents_.assign(num_vertices, 0);
    types_.assign(num_vertices, 0);
    // distances_with_vertices_.assign(num_vertices, std::make_pair(0,0));

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      std::swap(previous_order, current_order);
      InitializeDistances(sequence_graph, scoring_schema, sequence_bases[i],
                          previous_layer, previous_order, initialized_layer,
                          initialized_order);
      // InitializeDistancesWithSorting(sequence_bases[i], previous_layer,
      // previous_order, &initialized_layer, &initialized_order);
      current_layer = initialized_layer;
      PropagateInsertions(sequence_graph, scoring_schema, initialized_layer,
                          initialized_order, current_layer, current_order);
    }

    const ScoreType forward_alignment_cost = current_layer[current_order[0]];

    // For reverse complement.
    initialized_layer.assign(num_vertices, sequence_length * max_cost + 1);
    current_layer.assign(num_vertices, 0);

    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      current_order[i] = i;
    }

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      std::swap(previous_order, current_order);
      const char complement_base =
          sequence_batch.GetReverseComplementBaseOfSequenceAt(sequence_index,
                                                              i);
      InitializeDistances(sequence_graph, scoring_schema, complement_base,
                          previous_layer, previous_order, initialized_layer,
                          initialized_order);
      // InitializeDistancesWithSorting(base_complement_[sequence_bases[sequence_length
      // - 1 - i]], previous_layer, previous_order, &initialized_layer,
      // &initialized_order);
      current_layer = initialized_layer;
      PropagateInsertions(sequence_graph, scoring_schema, initialized_layer,
                          initialized_order, current_layer, current_order);
    }

    const ScoreType reverse_complement_alignment_cost =
        current_layer[current_order[0]];

    const ScoreType min_alignment_cost =
        std::min(forward_alignment_cost, reverse_complement_alignment_cost);

    std::cerr << "Sequence length: " << sequence_length
              << ", forward alignment cost:" << forward_alignment_cost
              << ", reverse complement alignment cost:"
              << reverse_complement_alignment_cost
              << ", alignment cost:" << min_alignment_cost << std::endl;
    return min_alignment_cost;
  }

 private:
  void PropagateInsertions(const SequenceGraph<GraphSizeType> &sequence_graph,
                           const ScoringSchema<ScoreType> &scoring_schema,
                           const std::vector<ScoreType> &initialized_layer,
                           const std::vector<GraphSizeType> &initialized_order,
                           std::vector<ScoreType> &current_layer,
                           std::vector<GraphSizeType> &current_order) {
    const GraphSizeType num_vertices = sequence_graph.GetNumVertices();
    GraphSizeType initialized_order_index = 0;
    GraphSizeType current_order_index = 0;
    visited_.assign(num_vertices, false);

    std::deque<GraphSizeType> updated_neighbors;

    while (initialized_order_index < num_vertices ||
           !updated_neighbors.empty()) {
      GraphSizeType min_vertex = num_vertices;

      if (initialized_order_index < num_vertices &&
          (updated_neighbors.empty() ||
           current_layer[initialized_order[initialized_order_index]] <
               current_layer[updated_neighbors.front()])) {
        min_vertex = initialized_order[initialized_order_index];
        ++initialized_order_index;
      } else {
        min_vertex = updated_neighbors.front();
        updated_neighbors.pop_front();
      }

      if (!visited_[min_vertex]) {
        visited_[min_vertex] = true;
        current_order[current_order_index] = min_vertex;
        ++current_order_index;
        for (const auto &neighbor : sequence_graph.GetNeighbors(min_vertex)) {
          if (!visited_[neighbor] &&
              current_layer[neighbor] > current_layer[min_vertex] +
                                            scoring_schema.insertion_penalty) {
            current_layer[neighbor] =
                current_layer[min_vertex] + scoring_schema.insertion_penalty;
            updated_neighbors.push_back(neighbor);
          }
        }
      }
    }
  }

  // Build the order look up table.
  void BuildOrderLookUpTable(const SequenceGraph<GraphSizeType> &sequence_graph,
                             const ScoringSchema<ScoreType> &scoring_schema,
                             const std::vector<ScoreType> &previous_layer,
                             const std::vector<GraphSizeType> &previous_order) {
    const GraphSizeType num_vertices = sequence_graph.GetNumVertices();
    GraphSizeType match_index = 0, substitution_index = 0, deletion_index = 0;
    GraphSizeType count = 0;

    while (count < 3 * num_vertices) {
      // Find the min.
      ScoreType min_distance =
          previous_layer[previous_order[num_vertices - 1]] +
          scoring_schema.substitution_penalty +
          scoring_schema.deletion_penalty + 1;
      int min_type = -1;  // 0 for match, 1 for substitution, 2 for deletion.
      if (match_index < num_vertices &&
          previous_layer[previous_order[match_index]] < min_distance) {
        min_distance = previous_layer[previous_order[match_index]];
        min_type = 0;
      }

      if (substitution_index < num_vertices &&
          previous_layer[previous_order[substitution_index]] +
                  scoring_schema.substitution_penalty <
              min_distance) {
        min_distance = previous_layer[previous_order[substitution_index]] +
                       scoring_schema.substitution_penalty;
        min_type = 1;
      }

      if (deletion_index < num_vertices &&
          previous_layer[previous_order[deletion_index]] +
                  scoring_schema.deletion_penalty <
              min_distance) {
        min_distance = previous_layer[previous_order[deletion_index]] +
                       scoring_schema.deletion_penalty;
        min_type = 2;
      }

      // Put the order of min into the look up table.
      if (min_type == 0) {
        order_look_up_table_[min_type * num_vertices +
                             previous_order[match_index]] = count;
        ++match_index;
      } else if (min_type == 1) {
        order_look_up_table_[min_type * num_vertices +
                             previous_order[substitution_index]] = count;
        ++substitution_index;
      } else {
        order_look_up_table_[min_type * num_vertices +
                             previous_order[deletion_index]] = count;
        ++deletion_index;
      }

      ++count;
    }
  }

  void InitializeDistancesWithSorting(
      const SequenceGraph<GraphSizeType> &sequence_graph,
      const ScoringSchema<ScoreType> &scoring_schema, const char sequence_base,
      const std::vector<ScoreType> &previous_layer,
      const std::vector<GraphSizeType> &previous_order,
      std::vector<ScoreType> &initialized_layer,
      std::vector<GraphSizeType> &initialized_order) {
    const GraphSizeType num_vertices = sequence_graph.GetNumVertices();

    // Initialize the layer.
    initialized_layer[0] = previous_layer[0] + scoring_schema.deletion_penalty;
    for (GraphSizeType j = 1; j < num_vertices; ++j) {
      ScoreType cost = 0;
      if (sequence_base != sequence_graph.GetVertexLabel(j)) {
        cost = scoring_schema.substitution_penalty;
      }
      initialized_layer[j] = previous_layer[0] + cost;
    }

    for (GraphSizeType i = 1; i < num_vertices; ++i) {
      if (initialized_layer[i] >
          previous_layer[i] + scoring_schema.deletion_penalty) {
        initialized_layer[i] =
            previous_layer[i] + scoring_schema.deletion_penalty;
      }

      for (const auto &neighbor : sequence_graph.GetNeighbors(i)) {
        ScoreType cost = 0;

        if (sequence_base != sequence_graph.GetVertexLabel(neighbor)) {
          cost = scoring_schema.substitution_penalty;
        }

        if (initialized_layer[neighbor] > previous_layer[i] + cost) {
          initialized_layer[neighbor] = previous_layer[i] + cost;
        }
      }
    }

    // Use sorting to get order.
    for (GraphSizeType vertex = 0; vertex < num_vertices; ++vertex) {
      distances_with_vertices_[vertex] =
          std::make_pair((*initialized_layer)[vertex], vertex);
    }
    std::sort(distances_with_vertices_.begin(), distances_with_vertices_.end());
    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      initialized_order[i] = distances_with_vertices_[i].second;
    }
  }

  void InitializeDistances(const SequenceGraph<GraphSizeType> &sequence_graph,
                           const ScoringSchema<ScoreType> &scoring_schema,
                           const char sequence_base,
                           const std::vector<ScoreType> &previous_layer,
                           const std::vector<GraphSizeType> &previous_order,
                           std::vector<ScoreType> &initialized_layer,
                           std::vector<GraphSizeType> &initialized_order) {
    const GraphSizeType num_vertices = sequence_graph.GetNumVertices();
    BuildOrderLookUpTable(sequence_graph, scoring_schema, previous_layer,
                          previous_order);

    // Initialize the layer
    initialized_layer[0] = previous_layer[0] + scoring_schema.deletion_penalty;
    parents_[0] = 0;
    types_[0] = 2;

    for (GraphSizeType j = 1; j < num_vertices; ++j) {
      ScoreType cost = 0;
      int type = 0;
      if (sequence_base != sequence_graph.GetVertexLabel(j)) {
        cost = scoring_schema.substitution_penalty;
        type = 1;
      }
      initialized_layer[j] = previous_layer[0] + cost;
      parents_[j] = 0;
      types_[j] = type;
    }

    for (GraphSizeType i = 1; i < num_vertices; ++i) {
      if (initialized_layer[i] >
          previous_layer[i] + scoring_schema.deletion_penalty) {
        initialized_layer[i] =
            previous_layer[i] + scoring_schema.deletion_penalty;
        parents_[i] = i;
        types_[i] = 2;
      }

      for (const auto &neighbor : sequence_graph.GetNeighbors(i)) {
        ScoreType cost = 0;
        int type = 0;

        if (sequence_base != sequence_graph.GetVertexLabel(neighbor)) {
          cost = scoring_schema.substitution_penalty;
          type = 1;
        }

        if (initialized_layer[neighbor] > previous_layer[i] + cost) {
          initialized_layer[neighbor] = previous_layer[i] + cost;
          parents_[neighbor] = i;
          types_[neighbor] = type;
        }
      }
    }

    // Get the order. One should notice there can be multiple vertices share the
    // same parent and type.
    order_offsets_.assign(3 * num_vertices + 1, 0);
    order_counts_.assign(3 * num_vertices, 0);

    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      order_offsets_
          [order_look_up_table_[types_[i] * num_vertices + parents_[i]] + 1]++;
    }

    for (GraphSizeType i = 1; i < 3 * num_vertices + 1; ++i) {
      order_offsets_[i] += order_offsets_[i - 1];
    }

    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      const GraphSizeType order =
          order_look_up_table_[types_[i] * num_vertices + parents_[i]];
      initialized_order[order_offsets_[order] + order_counts_[order]] = i;
      ++order_counts_[order];
    }
  }

  std::vector<GraphSizeType> order_look_up_table_;
  std::vector<bool> visited_;
  std::vector<GraphSizeType> parents_;
  std::vector<int> types_;
  std::vector<GraphSizeType> order_offsets_;
  std::vector<GraphSizeType> order_counts_;
  std::vector<std::pair<ScoreType, GraphSizeType>> distances_with_vertices_;
};

}  // namespace sgat
#endif  // SGAT_SGA_H
