#ifndef SGAT_DIJKSTRA_H
#define SGAT_DIJKSTRA_H

#include "scoring_schema.h"
#include "sequence_batch.h"
#include "sequence_graph.h"
#include "utils.h"
//#include "khash.h"

namespace sgat {

template <class GraphSizeType = int32_t, class QueryLengthType = int16_t,
          class ScoreType = int16_t>
struct VertexWithDistanceForDijkstra {
  GraphSizeType graph_vertex_id;
  QueryLengthType query_index;
  ScoreType distance;
  bool is_reverse_complementary;
};

template <class GraphSizeType = int32_t>
struct DijkstraAlgorithmStatistics {
  GraphSizeType forward_num_cells = 0;
  GraphSizeType rc_num_cells = 0;
};

// KHASH_MAP_INIT_INT(k32, uint32_t);

template <class GraphSizeType = int32_t, class QueryLengthType = int16_t,
          class ScoreType = int16_t>
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
            sequence_graph.GetNumVertices(), scoring_schema, stats);

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
  ScoreType AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
      uint32_t sequence_index, const SequenceBatch &sequence_batch,
      const SequenceGraph<GraphSizeType> &sequence_graph,
      GraphSizeType start_vertex,
      const ScoringSchema<ScoreType> &scoring_schema,
      DijkstraAlgorithmStatistics<GraphSizeType> &stats) {
    const GraphSizeType num_vertices = sequence_graph.GetNumVertices();
    const QueryLengthType sequence_length =
        sequence_batch.GetSequenceLengthAt(sequence_index);
    const std::string &sequence_bases =
        sequence_batch.GetSequenceAt(sequence_index);

    std::vector<std::unordered_map<QueryLengthType, ScoreType>>
        forward_vertex_distances(num_vertices);
    std::vector<std::unordered_map<QueryLengthType, ScoreType>>
        complementary_vertex_distances(num_vertices);

    // std::vector<khash_t(k32) *> forward_vertex_distances(num_vertices,
    // nullptr); std::vector<khash_t(k32) *>
    // complementary_vertex_distances(num_vertices,
    //                                                           nullptr);
    // for (GraphSizeType vi = 0; vi < num_vertices; ++vi) {
    //  forward_vertex_distances[vi] = kh_init(k32);
    //  complementary_vertex_distances[vi] = kh_init(k32);
    //}

    // std::vector<std::unordered_map<
    //    QueryLengthType, VertexWithDistanceForDijkstra<
    //                         GraphSizeType, QueryLengthType, ScoreType>>>
    //    vertex_parent(num_vertices);

    auto compare_function =
        [](const VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
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
                if (!v1.is_reverse_complementary &&
                    v2.is_reverse_complementary) {
                  return true;
                }
                return false;
              }
              return false;
            }
            return false;
          }
          return false;
        };

    std::priority_queue<VertexWithDistanceForDijkstra<
                            GraphSizeType, QueryLengthType, ScoreType>,
                        std::vector<VertexWithDistanceForDijkstra<
                            GraphSizeType, QueryLengthType, ScoreType>>,
                        decltype(compare_function)>
        Q(compare_function);

    stats.forward_num_cells = 0;
    stats.rc_num_cells = 0;

    GraphSizeType init_start_vertex =
        start_vertex == num_vertices ? 1 : start_vertex;
    GraphSizeType init_end_vertex =
        start_vertex == num_vertices ? num_vertices : start_vertex + 1;
    for (GraphSizeType vertex = init_start_vertex; vertex < init_end_vertex;
         ++vertex) {
      // Deal with forward strand fisrt.
      ScoreType cost = 0;
      const char vertex_label = sequence_graph.GetVertexLabel(vertex);

      if (sequence_bases[0] != vertex_label) {
        cost = scoring_schema.substitution_penalty;
      }
      cost = std::min(cost, scoring_schema.insertion_penalty);

      Q.push({/*graph_vertex_id=*/vertex, /*query_index=*/0, /*distance=*/cost,
              /*is_reverse_complementary=*/false});

      ++(stats.forward_num_cells);

      // int khash_return_code = 0;
      // khiter_t forward_vertex_distances_iterator =
      //    kh_put(k32, forward_vertex_distances[vertex], 0,
      //    &khash_return_code);
      // assert(khash_return_code != -1 && khash_return_code != 0);
      // kh_value(forward_vertex_distances[vertex],
      //         forward_vertex_distances_iterator) = cost;

      forward_vertex_distances[vertex][0] = cost;

      // Now deal with reverse complementary strand.
      cost = 0;
      const char complementary_vertex_label =
          sequence_graph.GetReverseComplementaryVertexLabel(vertex);

      if (sequence_bases[sequence_length - 1] != complementary_vertex_label) {
        cost = scoring_schema.substitution_penalty;
      }
      cost = std::min(cost, scoring_schema.insertion_penalty);

      Q.push({/*graph_vertex_id=*/vertex, /*query_index=*/0, /*distance=*/cost,
              /*is_reverse_complementary=*/true});

      // khash_return_code = 0;
      // khiter_t complementary_vertex_distances_iterator = kh_put(
      //    k32, complementary_vertex_distances[vertex], 0, &khash_return_code);
      // assert(khash_return_code != -1 && khash_return_code != 0);
      // kh_value(complementary_vertex_distances[vertex],
      //         complementary_vertex_distances_iterator) = cost;

      complementary_vertex_distances[vertex][0] = cost;
      ++(stats.rc_num_cells);

      // vertex_parent[vertex][0] =
      //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
      //                                  ScoreType>{
      //        /*graph_vertex_id=*/vertex, /*query_index=*/0,
      //        /*distance=*/cost};
      // std::cerr << "Init PUSH: " << vertex << " " << 0 << " " << cost
      //          << std::endl;
    }

    ScoreType min_alignment_cost = 0;

    while (!Q.empty()) {
      const auto current_vertex = Q.top();
      Q.pop();

      // Check if we reach the last layer where we can stop.
      if (current_vertex.query_index + 1 == sequence_length) {
        min_alignment_cost = current_vertex.distance;
        // auto previous_it =
        //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
        //                                  ScoreType>{-1, -1, 0};
        // auto it = current_vertex;
        // while (!(it.graph_vertex_id == previous_it.graph_vertex_id &&
        //         it.query_index == previous_it.query_index)) {
        //  std::cerr << "Traceback: "
        //            << "gi: " << it.graph_vertex_id << " qi: " <<
        //            it.query_index
        //            << " d: " << it.distance
        //            << " qb: " << sequence_bases[it.query_index]
        //            << " gb: " << labels_[it.graph_vertex_id];

        //  if (it.distance + substitution_penalty_ == previous_it.distance) {
        //    std::cerr << " op: M";
        //  }

        //  if (it.graph_vertex_id == previous_it.graph_vertex_id &&
        //      it.distance + deletion_penalty_ == previous_it.distance) {
        //    std::cerr << " op: D";
        //  }

        //  if (it.query_index == previous_it.query_index) {
        //    std::cerr << " op: I";
        //  }
        //  std::cerr << std::endl;
        //  previous_it = it;
        //  it = vertex_parent[it.graph_vertex_id][it.query_index];
        //}
        // std::cerr << std::endl;
        break;
      }

      auto &vertex_distances = current_vertex.is_reverse_complementary
                                   ? complementary_vertex_distances
                                   : forward_vertex_distances;

      const ScoreType min_insertion_distance =
          (ScoreType)(scoring_schema.insertion_penalty *
                      (current_vertex.query_index + 1));
      if (current_vertex.distance > min_insertion_distance) {
        Q.push(
            {/*graph_vertex_id=*/start_vertex,
             /*query_index=*/(QueryLengthType)(current_vertex.query_index + 1),
             /*distance=*/min_insertion_distance,
             /*is_reverse_complementary=*/
             current_vertex.is_reverse_complementary});
        vertex_distances[start_vertex][current_vertex.query_index + 1] =
            min_insertion_distance;
      }

      // Explore its neighbors.
      for (const auto &neighbor :
           sequence_graph.GetNeighbors(current_vertex.graph_vertex_id)) {
        // Process neighbors in the same layaer.
        const ScoreType new_deletion_distance =
            current_vertex.distance + scoring_schema.deletion_penalty;

        // khiter_t vertex_distances_iterator =
        //    kh_get(k32, vertex_distances[neighbor],
        //    current_vertex.query_index);
        // if (vertex_distances_iterator != kh_end(vertex_distances[neighbor]))
        // {
        //  if (new_deletion_distance <
        //      kh_value(vertex_distances[neighbor], vertex_distances_iterator))
        //      {
        //    kh_value(vertex_distances[neighbor], vertex_distances_iterator) =
        //        new_deletion_distance;
        if (vertex_distances[neighbor].find(current_vertex.query_index) !=
            vertex_distances[neighbor].end()) {
          if (new_deletion_distance <
              vertex_distances[neighbor][current_vertex.query_index]) {
            vertex_distances[neighbor][current_vertex.query_index] =
                new_deletion_distance;

            Q.push({/*graph_vertex_id=*/neighbor,
                    /*query_index=*/current_vertex.query_index,
                    /*distance=*/new_deletion_distance,
                    /*is_reverse_complementary=*/
                    current_vertex.is_reverse_complementary});
            // vertex_parent[neighbor][current_vertex.query_index] =
            //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
            //                                  ScoreType>{
            //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
            //        /*query_index=*/current_vertex.query_index,
            //        /*distance=*/current_vertex.distance};
            // std::cerr << "PUSH: " << neighbor << " "
            //          << current_vertex.query_index << " "
            //          << new_deletion_distance << std::endl;
            if (current_vertex.is_reverse_complementary) {
              ++(stats.rc_num_cells);
            } else {
              ++(stats.forward_num_cells);
            }
          }
        } else {
          // int khash_return_code = 0;
          // khiter_t vertex_distances_iterator =
          //    kh_put(k32, vertex_distances[neighbor],
          //           current_vertex.query_index, &khash_return_code);
          // assert(khash_return_code != -1 && khash_return_code != 0);
          // kh_value(vertex_distances[neighbor], vertex_distances_iterator) =
          //    new_deletion_distance;

          vertex_distances[neighbor][current_vertex.query_index] =
              new_deletion_distance;
          Q.push({/*graph_vertex_id=*/neighbor,
                  /*query_index=*/current_vertex.query_index,
                  /*distance=*/new_deletion_distance,
                  /*is_reverse_complementary=*/
                  current_vertex.is_reverse_complementary});
          // vertex_parent[neighbor][current_vertex.query_index] =
          //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
          //                                  ScoreType>{
          //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
          //        /*query_index=*/current_vertex.query_index,
          //        /*distance=*/current_vertex.distance};

          // std::cerr << "PUSH: " << neighbor << " " <<
          // current_vertex.query_index
          //          << " " << new_deletion_distance << std::endl;
          if (current_vertex.is_reverse_complementary) {
            ++(stats.rc_num_cells);
          } else {
            ++(stats.forward_num_cells);
          }
        }

        // Process neighbors in the next layaer.
        const QueryLengthType query_index = current_vertex.query_index + 1;

        ScoreType cost = 0;
        const char vertex_label =
            current_vertex.is_reverse_complementary
                ? sequence_graph.GetReverseComplementaryVertexLabel(neighbor)
                : sequence_graph.GetVertexLabel(neighbor);
        const char sequence_base =
            current_vertex.is_reverse_complementary
                ? sequence_bases[sequence_length - 1 - query_index]
                : sequence_bases[query_index];

        if (sequence_base != vertex_label) {
          cost = scoring_schema.substitution_penalty;
        }

        const ScoreType new_match_or_mismatch_distance =
            std::min(
                current_vertex.distance,
                (ScoreType)(scoring_schema.insertion_penalty * query_index)) +
            cost;

        // std::cerr << "seq base: " << sequence_bases[query_index] << " label:
        // " << labels_[neighbor] << " cost: " << cost << " d: " <<
        // new_match_or_mismatch_distance << std::endl;

        // vertex_distances_iterator =
        //    kh_get(k32, vertex_distances[neighbor], query_index);
        // if (vertex_distances_iterator != kh_end(vertex_distances[neighbor]))
        // {
        //  if (new_match_or_mismatch_distance <
        //      kh_value(vertex_distances[neighbor], vertex_distances_iterator))
        //      {
        //    kh_value(vertex_distances[neighbor], vertex_distances_iterator) =
        //        new_match_or_mismatch_distance;

        if (vertex_distances[neighbor].find(query_index) !=
            vertex_distances[neighbor].end()) {
          if (new_match_or_mismatch_distance <
              vertex_distances[neighbor][query_index]) {
            vertex_distances[neighbor][query_index] =
                new_match_or_mismatch_distance;

            Q.push({/*graph_vertex_id=*/neighbor,
                    /*query_index=*/query_index,
                    /*distance=*/new_match_or_mismatch_distance,
                    /*is_reverse_complementary=*/
                    current_vertex.is_reverse_complementary});
            // vertex_parent[neighbor][query_index] =
            //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
            //                                  ScoreType>{
            //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
            //        /*query_index=*/current_vertex.query_index,
            //        /*distance=*/current_vertex.distance};

            // std::cerr << "PUSH: " << neighbor << " " << query_index << " "
            //          << new_match_or_mismatch_distance << std::endl;

            if (current_vertex.is_reverse_complementary) {
              ++(stats.rc_num_cells);
            } else {
              ++(stats.forward_num_cells);
            }
          }
        } else {
          // int khash_return_code = 0;
          // khiter_t vertex_distances_iterator = kh_put(
          //    k32, vertex_distances[neighbor], query_index,
          //    &khash_return_code);
          // assert(khash_return_code != -1 && khash_return_code != 0);
          // kh_value(vertex_distances[neighbor], vertex_distances_iterator) =
          //    new_match_or_mismatch_distance;

          vertex_distances[neighbor][query_index] =
              new_match_or_mismatch_distance;

          Q.push({/*graph_vertex_id=*/neighbor,
                  /*query_index=*/query_index,
                  /*distance=*/new_match_or_mismatch_distance,
                  /*is_reverse_complementary=*/
                  current_vertex.is_reverse_complementary});
          // vertex_parent[neighbor][query_index] =
          //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
          //                                  ScoreType>{
          //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
          //        /*query_index=*/current_vertex.query_index,
          //        /*distance=*/current_vertex.distance};

          // std::cerr << "PUSH: " << neighbor << " " << query_index << " "
          //          << new_match_or_mismatch_distance << std::endl;
          if (current_vertex.is_reverse_complementary) {
            ++(stats.rc_num_cells);
          } else {
            ++(stats.forward_num_cells);
          }
        }
      }  // End for exploring neighbor.

      // Process insertions to the next layaer.
      const QueryLengthType query_index = current_vertex.query_index + 1;

      const ScoreType new_insertion_distance =
          current_vertex.distance + scoring_schema.insertion_penalty;

      // khiter_t vertex_distances_iterator = kh_get(
      //    k32, vertex_distances[current_vertex.graph_vertex_id], query_index);
      // if (vertex_distances_iterator !=
      //    kh_end(vertex_distances[current_vertex.graph_vertex_id])) {
      //  if (new_insertion_distance <
      //      kh_value(vertex_distances[current_vertex.graph_vertex_id],
      //               vertex_distances_iterator)) {
      //    kh_value(vertex_distances[current_vertex.graph_vertex_id],
      //             vertex_distances_iterator) = new_insertion_distance;

      if (vertex_distances[current_vertex.graph_vertex_id].find(query_index) !=
          vertex_distances[current_vertex.graph_vertex_id].end()) {
        if (new_insertion_distance <
            vertex_distances[current_vertex.graph_vertex_id][query_index]) {
          vertex_distances[current_vertex.graph_vertex_id][query_index] =
              new_insertion_distance;

          Q.push({/*graph_vertex_id=*/current_vertex.graph_vertex_id,
                  /*query_index=*/query_index,
                  /*distance=*/new_insertion_distance,
                  /*is_reverse_complementary=*/
                  current_vertex.is_reverse_complementary});
          // vertex_parent[current_vertex.graph_vertex_id][query_index] =
          //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
          //                                  ScoreType>{
          //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
          //        /*query_index=*/current_vertex.query_index,
          //        /*distance=*/current_vertex.distance};

          // std::cerr << "PUSH: " << current_vertex.graph_vertex_id << " "
          //          << query_index << " " << new_insertion_distance
          //          << std::endl;
          if (current_vertex.is_reverse_complementary) {
            ++(stats.rc_num_cells);
          } else {
            ++(stats.forward_num_cells);
          }
        }
      } else {
        // int khash_return_code = 0;
        // khiter_t vertex_distances_iterator =
        //    kh_put(k32, vertex_distances[current_vertex.graph_vertex_id],
        //           query_index, &khash_return_code);
        // assert(khash_return_code != -1 && khash_return_code != 0);
        // kh_value(vertex_distances[current_vertex.graph_vertex_id],
        //         vertex_distances_iterator) = new_insertion_distance;

        vertex_distances[current_vertex.graph_vertex_id][query_index] =
            new_insertion_distance;

        Q.push({/*graph_vertex_id=*/current_vertex.graph_vertex_id,
                /*query_index=*/query_index,
                /*distance=*/new_insertion_distance,
                /*is_reverse_complementary=*/
                current_vertex.is_reverse_complementary});
        // vertex_parent[current_vertex.graph_vertex_id][query_index] =
        //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
        //                                  ScoreType>{
        //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
        //        /*query_index=*/current_vertex.query_index,
        //        /*distance=*/current_vertex.distance};

        // std::cerr << "PUSH: " << current_vertex.graph_vertex_id << " "
        //          << query_index << " " << new_insertion_distance <<
        //          std::endl;
        if (current_vertex.is_reverse_complementary) {
          ++(stats.rc_num_cells);
        } else {
          ++(stats.forward_num_cells);
        }
      }
    }

    // for (GraphSizeType vi = 0; vi < num_vertices; ++vi) {
    //  kh_destroy(k32, forward_vertex_distances[vi]);
    //  kh_destroy(k32, complementary_vertex_distances[vi]);
    //}

    // std::cerr << "Sequence length: " << sequence_length
    //          << ", forward alignment cost:" << forward_alignment_cost
    //          << ", reverse complement alignment cost:"
    //          << reverse_complement_alignment_cost
    //          << ", alignment cost:" << min_alignment_cost << std::endl;
    return min_alignment_cost;
  }
};

}  // namespace sgat
#endif  // SGAT_DIJKSTRA_H
