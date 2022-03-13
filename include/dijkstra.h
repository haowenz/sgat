#ifndef SGAT_DIJKSTRA_H
#define SGAT_DIJKSTRA_H

#include "scoring_schema.h"
#include "sequence_batch.h"
#include "sequence_graph.h"
#include "utils.h"

namespace sgat {

template <class GraphSizeType, class QueryLengthType, class ScoreType>
struct VertexWithDistanceForDijkstra {
  GraphSizeType graph_vertex_id;
  GraphSizeType graph_vertex_j;
  QueryLengthType query_index;
  ScoreType distance;
  bool is_reverse_complementary;
};

template <class GraphSizeType, class QueryLengthType, class ScoreType>
bool CompareVertexWithDistanceForDijkstra(
    const VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
                                        ScoreType> &v1,
    const VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
                                        ScoreType> &v2) {
  if (v1.distance > v2.distance) {
    return true;
  }

  if (v1.distance == v2.distance) {
    if (v1.query_index < v2.query_index) {
      return true;
    }
    if (v1.query_index == v2.query_index) {
      if (v1.graph_vertex_id > v2.graph_vertex_id) {
        return true;
      }
      if (v1.graph_vertex_id == v2.graph_vertex_id) {
        if (v1.graph_vertex_j > v2.graph_vertex_j) {
          return true;
        }
        if (v1.graph_vertex_j == v2.graph_vertex_j) {
          if (!v1.is_reverse_complementary && v2.is_reverse_complementary) {
            return true;
          }
          return false;
        }
        return false;
      }
      return false;
    }
    return false;
  }
  return false;
}

template <class GraphSizeType>
struct DijkstraAlgorithmStatistics {
  GraphSizeType forward_num_cells = 0;
  GraphSizeType rc_num_cells = 0;
};

template <class GraphSizeType, class QueryLengthType, class ScoreType>
class DijkstraAligner {
 public:
  DijkstraAligner() = default;
  ~DijkstraAligner() = default;

  ScoreType AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
      uint32_t sequence_index, const SequenceBatch &sequence_batch,
      const SequenceGraph<GraphSizeType> &sequence_graph,
      const ScoringSchema<ScoreType> &scoring_schema) {
    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    DijkstraAlgorithmStatistics<GraphSizeType> stats;
    const ScoreType min_alignment_cost =
        AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
            sequence_index, sequence_batch, sequence_graph,
            sequence_graph.GetNumVerticesInCompactedGraph(), scoring_schema,
            stats);

    std::cerr << "Sequence length: " << sequence_length
              << ", alignment cost:" << min_alignment_cost
              << ", forward num cells:" << stats.forward_num_cells
              << ", reverse num cells: " << stats.rc_num_cells << std::endl;
    return min_alignment_cost;
  }

  ScoreType ExtendUsingLinearGapPenaltyWithDijkstraAlgorithm(
      uint32_t sequence_index, const SequenceBatch &sequence_batch,
      const SequenceGraph<GraphSizeType> &sequence_graph,

      GraphSizeType start_vertex,
      const ScoringSchema<ScoreType> &scoring_schema) {
    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    DijkstraAlgorithmStatistics<GraphSizeType> stats;
    const ScoreType min_alignment_cost =
        AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
            sequence_index, sequence_batch, sequence_graph, start_vertex,
            scoring_schema, stats);

    std::cerr << "Sequence length: " << sequence_length
              << ", alignment cost:" << min_alignment_cost
              << ", forward num cells:" << stats.forward_num_cells
              << ", reverse num cells: " << stats.rc_num_cells << std::endl;
    return min_alignment_cost;
  }

 private:
  // Return true if it got updated.
  bool UpdateVertexDistanceIfNecessary(
      GraphSizeType vertex_id, GraphSizeType vertex_j,
      QueryLengthType query_index, ScoreType distance,
      std::vector<std::vector<std::unordered_map<QueryLengthType, ScoreType>>>
          &vertex_distances) {
    const bool is_entry_found =
        vertex_distances[vertex_id][vertex_j].find(query_index) !=
        vertex_distances[vertex_id][vertex_j].end();

    bool is_updated = false;
    if ((is_entry_found &&
         distance < vertex_distances[vertex_id][vertex_j][query_index]) ||
        !is_entry_found) {
      vertex_distances[vertex_id][vertex_j][query_index] = distance;
      is_updated = true;
    }

    return is_updated;
  }

  // Return true if it got updated.
  bool PushDijkstraQueueAndUpdateDistanceIfNecessary(
      GraphSizeType vertex_id, GraphSizeType vertex_j,
      QueryLengthType query_index, ScoreType vertex_distance,
      bool is_reverse_complementary,
      std::vector<std::vector<std::unordered_map<QueryLengthType, ScoreType>>>
          &vertex_distances,
      std::priority_queue<VertexWithDistanceForDijkstra<
                              GraphSizeType, QueryLengthType, ScoreType>,
                          std::vector<VertexWithDistanceForDijkstra<
                              GraphSizeType, QueryLengthType, ScoreType>>,
                          decltype(&CompareVertexWithDistanceForDijkstra<
                                   GraphSizeType, QueryLengthType, ScoreType>)>
          &Q,
      DijkstraAlgorithmStatistics<GraphSizeType> &stats) {
    const bool is_updated = UpdateVertexDistanceIfNecessary(
        vertex_id, vertex_j, query_index, vertex_distance, vertex_distances);
    if (is_updated) {
      Q.push({/*graph_vertex_id=*/vertex_id,
              /*graph_vertex_j=*/vertex_j,
              /*query_index=*/query_index,
              /*distance=*/vertex_distance,
              /*is_reverse_complementary=*/
              is_reverse_complementary});

      if (is_reverse_complementary) {
        ++(stats.rc_num_cells);
      } else {
        ++(stats.forward_num_cells);
      }
    }

    return is_updated;
  }

  void ProcessNeighbors(
      GraphSizeType neighbor,
      const VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
                                          ScoreType> &current_vertex,
      uint32_t sequence_index, const SequenceBatch &sequence_batch,
      GraphSizeType start_vertex,
      const SequenceGraph<GraphSizeType> &sequence_graph,
      const ScoringSchema<ScoreType> &scoring_schema,
      std::vector<std::vector<std::unordered_map<QueryLengthType, ScoreType>>>
          &vertex_distances,
      std::priority_queue<VertexWithDistanceForDijkstra<
                              GraphSizeType, QueryLengthType, ScoreType>,
                          std::vector<VertexWithDistanceForDijkstra<
                              GraphSizeType, QueryLengthType, ScoreType>>,
                          decltype(&CompareVertexWithDistanceForDijkstra<
                                   GraphSizeType, QueryLengthType, ScoreType>)>
          &Q,
      DijkstraAlgorithmStatistics<GraphSizeType> &stats) {
    const GraphSizeType num_vertices =
        sequence_graph.GetNumVerticesInCompactedGraph();

    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    const std::string &sequence_bases =
        sequence_batch.GetSequenceAt(sequence_index);

    // This is the last char in the vertex, we have to add its neighbors to
    // the queue.
    // Add neighbors in the same layaer.
    const ScoreType new_deletion_distance =
        current_vertex.distance + scoring_schema.deletion_penalty;

    PushDijkstraQueueAndUpdateDistanceIfNecessary(
        neighbor, /*vertex_j=*/0, current_vertex.query_index,
        new_deletion_distance, current_vertex.is_reverse_complementary,
        vertex_distances, Q, stats);

    // Add neighbors in the next layaer.
    const char vertex_0_label =
        current_vertex.is_reverse_complementary
            ? sequence_graph.GetComplementaryVertexLabelInCompactedGraph(
                  neighbor, 0)
            : sequence_graph.GetVertexLabelInCompactedGraph(neighbor, 0);
    const char sequence_base =
        current_vertex.is_reverse_complementary
            ? sequence_bases[sequence_length - 1 -
                             (current_vertex.query_index + 1)]
            : sequence_bases[(current_vertex.query_index + 1)];

    const bool is_mismatch = (vertex_0_label != sequence_base);

    const ScoreType new_substitution_distance =
        current_vertex.distance +
        (is_mismatch ? scoring_schema.substitution_penalty : 0);

    PushDijkstraQueueAndUpdateDistanceIfNecessary(
        neighbor, 0, current_vertex.query_index + 1, new_substitution_distance,
        current_vertex.is_reverse_complementary, vertex_distances, Q, stats);
  }

  ScoreType AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
      uint32_t sequence_index, const SequenceBatch &sequence_batch,
      const SequenceGraph<GraphSizeType> &sequence_graph,
      GraphSizeType start_vertex,
      const ScoringSchema<ScoreType> &scoring_schema,
      DijkstraAlgorithmStatistics<GraphSizeType> &stats) {
    const GraphSizeType num_vertices =
        sequence_graph.GetNumVerticesInCompactedGraph();
    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    const std::string &sequence_bases =
        sequence_batch.GetSequenceAt(sequence_index);

    std::vector<std::vector<std::unordered_map<QueryLengthType, ScoreType>>>
        forward_vertex_distances;
    std::vector<std::vector<std::unordered_map<QueryLengthType, ScoreType>>>
        complementary_vertex_distances;

    for (GraphSizeType vertex_id = 0; vertex_id < num_vertices; ++vertex_id) {
      GraphSizeType vertex_length =
          sequence_graph.GetVertexLengthInCompactedGraph(vertex_id);
      forward_vertex_distances.emplace_back(
          std::vector<std::unordered_map<QueryLengthType, ScoreType>>(
              vertex_length));
      complementary_vertex_distances.emplace_back(
          std::vector<std::unordered_map<QueryLengthType, ScoreType>>(
              vertex_length));
    }

    std::priority_queue<VertexWithDistanceForDijkstra<
                            GraphSizeType, QueryLengthType, ScoreType>,
                        std::vector<VertexWithDistanceForDijkstra<
                            GraphSizeType, QueryLengthType, ScoreType>>,
                        decltype(&CompareVertexWithDistanceForDijkstra<
                                 GraphSizeType, QueryLengthType, ScoreType>)>
        Q(CompareVertexWithDistanceForDijkstra);

    stats.forward_num_cells = 0;
    stats.rc_num_cells = 0;

    // Initialization.
    PushDijkstraQueueAndUpdateDistanceIfNecessary(
        0, 0, -1, /*distance=*/0, false, forward_vertex_distances, Q, stats);

    if (start_vertex == num_vertices) {
      PushDijkstraQueueAndUpdateDistanceIfNecessary(
          0, 0, -1, /*distance=*/0, true, complementary_vertex_distances, Q,
          stats);

      for (GraphSizeType vertex_id = 1; vertex_id < num_vertices; ++vertex_id) {
        const GraphSizeType vertex_length =
            sequence_graph.GetVertexLengthInCompactedGraph(vertex_id);
        for (GraphSizeType vertex_j = 0; vertex_j < vertex_length; ++vertex_j) {
          PushDijkstraQueueAndUpdateDistanceIfNecessary(
              vertex_id, vertex_j, -1, /*distance=*/0, false,
              forward_vertex_distances, Q, stats);

          PushDijkstraQueueAndUpdateDistanceIfNecessary(
              vertex_id, vertex_j, -1, /*distance=*/0, true,
              complementary_vertex_distances, Q, stats);
        }
      }
    }

    ScoreType min_alignment_cost = 0;

    while (!Q.empty()) {
      const auto current_vertex = Q.top();
      Q.pop();

#ifdef SGAT_DEBUG
      std::cerr << "vi: " << current_vertex.graph_vertex_id
                << ", vj: " << current_vertex.graph_vertex_j
                << ", qi: " << current_vertex.query_index
                << ", d: " << current_vertex.distance
                << ", rc: " << current_vertex.is_reverse_complementary
                << std::endl;
#endif

      // Check if we reach the last layer where we can stop.
      if (current_vertex.query_index + 1 == sequence_length) {
        min_alignment_cost = current_vertex.distance;
        break;
      }

      const GraphSizeType vertex_id = current_vertex.graph_vertex_id;
      const GraphSizeType vertex_j = current_vertex.graph_vertex_j;
      const QueryLengthType query_index = current_vertex.query_index;
      const ScoreType vertex_distance = current_vertex.distance;
      const bool is_reverse_complementary =
          current_vertex.is_reverse_complementary;
      const GraphSizeType vertex_length =
          sequence_graph.GetVertexLengthInCompactedGraph(vertex_id);

      auto &vertex_distances = is_reverse_complementary
                                   ? complementary_vertex_distances
                                   : forward_vertex_distances;

      // Process insertions in next layer.
      const ScoreType new_insertion_distance =
          vertex_distance + scoring_schema.insertion_penalty;

      PushDijkstraQueueAndUpdateDistanceIfNecessary(
          vertex_id, vertex_j, query_index + 1, new_insertion_distance,
          is_reverse_complementary, vertex_distances, Q, stats);

      // Deal with the case when the position is not at the end of a vertex.
      if (vertex_j < vertex_length - 1) {
        const ScoreType new_deletion_distance =
            vertex_distance + scoring_schema.deletion_penalty;
        PushDijkstraQueueAndUpdateDistanceIfNecessary(
            vertex_id, vertex_j + 1, query_index, new_deletion_distance,
            is_reverse_complementary, vertex_distances, Q, stats);

        const char vertex_next_label =
            current_vertex.is_reverse_complementary
                ? sequence_graph.GetComplementaryVertexLabelInCompactedGraph(
                      vertex_id, vertex_j + 1)
                : sequence_graph.GetVertexLabelInCompactedGraph(vertex_id,
                                                                vertex_j + 1);
        const char sequence_base =
            current_vertex.is_reverse_complementary
                ? sequence_bases[sequence_length - 1 - (query_index + 1)]
                : sequence_bases[query_index + 1];

        const bool is_mismatch = (vertex_next_label != sequence_base);

        // We can deal with the match immediately if we want.
        const ScoreType new_substitution_distance =
            vertex_distance +
            (is_mismatch ? scoring_schema.substitution_penalty : 0);

        PushDijkstraQueueAndUpdateDistanceIfNecessary(
            vertex_id, vertex_j + 1, query_index + 1, new_substitution_distance,
            is_reverse_complementary, vertex_distances, Q, stats);

        continue;
      }

      // This is the last char in the vertex, we have to add its neighbors to
      // the queue.
      if (vertex_id == 0) {
        const GraphSizeType start_neighbor =
            start_vertex == num_vertices ? 1 : start_vertex;
        const GraphSizeType end_neighbor =
            start_vertex == num_vertices ? num_vertices : start_vertex + 1;
        for (GraphSizeType neighbor = start_neighbor; neighbor < end_neighbor;
             ++neighbor) {
          ProcessNeighbors(neighbor, current_vertex, sequence_index,
                           sequence_batch, start_vertex, sequence_graph,
                           scoring_schema, vertex_distances, Q, stats);
        }
        continue;
      }

      for (const auto &neighbor :
           sequence_graph.GetNeighborsInCompatedGraph(vertex_id)) {
        ProcessNeighbors(neighbor, current_vertex, sequence_index,
                         sequence_batch, start_vertex, sequence_graph,
                         scoring_schema, vertex_distances, Q, stats);
      }
    }  // End while.

    return min_alignment_cost;
  }
};

}  // namespace sgat
#endif  // SGAT_DIJKSTRA_H
