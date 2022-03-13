#ifndef SGAT_SEQUENCEGRAPH_H
#define SGAT_SEQUENCEGRAPH_H

#include <assert.h>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "gfa.h"
#include "sequence_batch.h"
#include "utils.h"

#define SGAT_VERSION "0.1.0-r3"

namespace sgat {

template <class GraphSizeType = int32_t>
class SequenceGraph {
 public:
  SequenceGraph() = default;
  ~SequenceGraph() = default;

  inline GraphSizeType GetNumVerticesInCompactedGraph() const {
    return compacted_graph_labels_.size();
  }

  inline GraphSizeType GetVertexLengthInCompactedGraph(
      GraphSizeType vertex_id) const {
    return compacted_graph_labels_[vertex_id].length();
  }

  inline char GetVertexLabelInCompactedGraph(GraphSizeType vertex_id,
                                             GraphSizeType j) const {
    return compacted_graph_labels_[vertex_id][j];
  }

  inline GraphSizeType GetNumEdgesInCompactedGraph() const {
    GraphSizeType num_edges = 0;
    for (const std::vector<GraphSizeType> &neighbors :
         compacted_graph_adjacency_list_) {
      num_edges += neighbors.size();
    }
    return num_edges;
  }

  inline const GraphSizeType GetNumVertices() const { return labels_.size(); }

  GraphSizeType GetNumEdges() const {
    GraphSizeType num_edges = 0;
    for (const std::vector<GraphSizeType> &neighbors : adjacency_list_) {
      num_edges += neighbors.size();
    }
    return num_edges;
  }

  void AddReverseComplementaryVertexIfNecessary(
      const gfa_t *gfa_graph, uint32_t gfa_vertex_id,
      std::vector<GraphSizeType> &reverse_complementary_compacted_vertex_id) {
    const uint32_t gfa_segment_id = gfa_vertex_id >> 1;
    const GraphSizeType vertex_id = gfa_segment_id + 1;
    const bool is_gfa_segment_forward = ((gfa_vertex_id & 1) == 0);
    const bool vertex_rc_is_not_added =
        (reverse_complementary_compacted_vertex_id[vertex_id] == vertex_id);

    // std::cerr << "gfa_vi: " << gfa_vertex_id << " si: " << gfa_segment_id
    //          << " rcvi: "
    //          << reverse_complementary_compacted_vertex_id[vertex_id]
    //          << " vi: " << vertex_id << " forward:" << is_gfa_segment_forward
    //          << std::endl;
    if (!is_gfa_segment_forward && vertex_rc_is_not_added) {
      const gfa_seg_t &gfa_segment = gfa_graph->seg[gfa_segment_id];
      std::string rc_sequence_bases;
      for (GraphSizeType qi = 0; qi < (GraphSizeType)gfa_segment.len; ++qi) {
        rc_sequence_bases.push_back(
            base_complement_[(int)gfa_segment.seq
                                 [(GraphSizeType)gfa_segment.len - 1 - qi]]);
      }

      compacted_graph_labels_.emplace_back(rc_sequence_bases);
      compacted_graph_adjacency_list_.emplace_back(
          std::vector<GraphSizeType>());

      const uint32_t new_compacted_graph_rc_vertex_id =
          compacted_graph_adjacency_list_.size() - 1;
      reverse_complementary_compacted_vertex_id[vertex_id] =
          new_compacted_graph_rc_vertex_id;
    }
  }

  void AddEdge(const gfa_t *gfa_graph, const gfa_arc_t &gfa_arc,
               const std::vector<GraphSizeType>
                   &reverse_complementary_compacted_vertex_id) {
    const uint32_t gfa_head_vertex_id = gfa_arc_head(gfa_arc);
    const uint32_t gfa_head_segment_id = gfa_head_vertex_id >> 1;
    const bool is_gfa_head_segment_forward = ((gfa_head_vertex_id & 1) == 0);

    GraphSizeType head_vertex_id = gfa_head_segment_id + 1;
    if (!is_gfa_head_segment_forward) {
      head_vertex_id =
          reverse_complementary_compacted_vertex_id[head_vertex_id];
    }

    const uint32_t gfa_tail_vertex_id = gfa_arc_tail(gfa_arc);
    const uint32_t gfa_tail_segment_id = gfa_tail_vertex_id >> 1;
    const bool is_gfa_tail_segment_forward = ((gfa_tail_vertex_id & 1) == 0);

    GraphSizeType tail_vertex_id = gfa_tail_segment_id + 1;
    if (!is_gfa_tail_segment_forward) {
      tail_vertex_id =
          reverse_complementary_compacted_vertex_id[tail_vertex_id];
    }

    compacted_graph_adjacency_list_[head_vertex_id].emplace_back(
        tail_vertex_id);
  }

  // An algorithm to convert a GFA to a compacted sequence graph with two passes
  // on arcs and one pass on segments. We assume there is no overlap. The
  // algorithm first passes the arcs to know which node needs to be duplicated
  // for a reverse complement. Then the segments are parsed to create the
  // vertices and labels. Finally the arcs are parsed again to add the edges.
  void LoadFromGfaFile(const std::string &graph_file_path) {
    gfa_t *gfa_graph = gfa_read(graph_file_path.data());
    // gfa_print(gfa_graph, stderr, 0);

    std::vector<GraphSizeType> reverse_complementary_compacted_vertex_id;

    const uint32_t num_segments = gfa_graph->n_seg;
    reverse_complementary_compacted_vertex_id.reserve(num_segments);
    compacted_graph_labels_.reserve(num_segments);
    compacted_graph_adjacency_list_.reserve(num_segments);

    // Add a dummy vertex.
    compacted_graph_labels_.emplace_back("N");
    reverse_complementary_compacted_vertex_id.emplace_back(0);
    compacted_graph_adjacency_list_.emplace_back(std::vector<GraphSizeType>());

    for (uint32_t si = 0; si < num_segments; ++si) {
      const gfa_seg_t &current_segment = gfa_graph->seg[si];
      compacted_graph_labels_.emplace_back(std::string(current_segment.seq));

      const uint32_t compacted_graph_vertex_id = si + 1;
      reverse_complementary_compacted_vertex_id.emplace_back(
          compacted_graph_vertex_id);

      compacted_graph_adjacency_list_.emplace_back(
          std::vector<GraphSizeType>());
    }

    const uint64_t num_arcs = gfa_graph->n_arc;
    // std::cerr << "NUM ARCS: " << num_arcs << std::endl;
    const gfa_arc_t *gfa_arcs = gfa_graph->arc;

    for (uint64_t ai = 0; ai < num_arcs; ++ai) {
      if (gfa_arcs[ai].del || gfa_arcs[ai].comp) continue;
      const uint32_t gfa_head_vertex_id = gfa_arc_head(gfa_arcs[ai]);
      AddReverseComplementaryVertexIfNecessary(
          gfa_graph, gfa_head_vertex_id,
          reverse_complementary_compacted_vertex_id);

      const uint32_t gfa_tail_vertex_id = gfa_arc_tail(gfa_arcs[ai]);
      AddReverseComplementaryVertexIfNecessary(
          gfa_graph, gfa_tail_vertex_id,
          reverse_complementary_compacted_vertex_id);
    }

    for (uint64_t ai = 0; ai < num_arcs; ++ai) {
      if (gfa_arcs[ai].del || gfa_arcs[ai].comp) continue;
      AddEdge(gfa_graph, gfa_arcs[ai],
              reverse_complementary_compacted_vertex_id);
    }
  }

  void LoadFromTxtFile(const std::string &graph_file_path) {
    std::string line;
    std::ifstream infile(graph_file_path);

    GraphSizeType num_vertices = 0;
    GraphSizeType row_index = 0;

    compacted_graph_labels_.emplace_back("N");
    compacted_graph_adjacency_list_.emplace_back(std::vector<GraphSizeType>());

    while (std::getline(infile, line)) {
      std::istringstream inputString(line);
      // get count of vertices from header row
      if (row_index == 0) {
        inputString >> num_vertices;
        compacted_graph_labels_.reserve(num_vertices);
      } else {  // get out-neighbor vertex ids and vertex label
        assert(row_index <= num_vertices);
        compacted_graph_adjacency_list_.emplace_back(
            std::vector<GraphSizeType>());

        // Parse the input line
        std::vector<std::string> tokens(
            std::istream_iterator<std::string>{inputString},
            std::istream_iterator<std::string>());
        assert(tokens.size() > 0);
        compacted_graph_labels_.emplace_back(tokens.back());
        for (auto it = tokens.begin();
             it != tokens.end() && std::next(it) != tokens.end(); it++) {
          compacted_graph_adjacency_list_.back().emplace_back(stoi(*it) + 1);
        }
      }
      row_index++;
    }

    std::cerr << "# vertices in compacted graph: "
              << GetNumVerticesInCompactedGraph()
              << ", # edges in compacted graph: "
              << GetNumEdgesInCompactedGraph() << std::endl;
  }

  void OutputCompactedGraphInGFA(std::string &output_file_path) const {
    std::ofstream outstrm(output_file_path);
    outstrm << "H\tVN:Z:1.0\n";
    // Both for loops start from 1 to skip the dummy vertex.
    for (uint32_t i = 1; i < compacted_graph_labels_.size(); ++i) {
      outstrm << "S\t" << i << "\t" << compacted_graph_labels_[i] << "\n";
    }

    for (uint32_t i = 1; i < compacted_graph_adjacency_list_.size(); ++i) {
      for (auto neighbor : compacted_graph_adjacency_list_[i]) {
        outstrm << "L\t" << i << "\t+\t" << neighbor << "\t+\t0M\n";
      }
    }
  }

  void OutputCharLabeledGraphInGFA(std::string &output_file_path) const {
    std::ofstream outstrm(output_file_path);
    outstrm << "H\tVN:Z:1.0\n";
    // Both for loops start from 1 to skip the dummy vertex.
    for (uint32_t i = 1; i < labels_.size(); ++i) {
      outstrm << "S\t" << i << "\t" << labels_[i] << "\n";
    }

    for (uint32_t i = 1; i < adjacency_list_.size(); ++i) {
      for (auto neighbor : adjacency_list_[i]) {
        outstrm << "L\t" << i << "\t+\t" << neighbor << "\t+\t0M\n";
      }
    }
  }

  void GenerateCharLabeledGraph() {
    for (const std::string &compacted_graph_label : compacted_graph_labels_) {
      labels_.emplace_back(compacted_graph_label[0]);
      adjacency_list_.emplace_back(std::vector<GraphSizeType>());
    }

    // Keep the original vertex ids unchanged.
    GraphSizeType vertex_id = compacted_graph_labels_.size();
    GraphSizeType compacted_graph_vertex_id = 0;

    for (const std::string &compacted_graph_label : compacted_graph_labels_) {
      GraphSizeType compacted_graph_label_length =
          compacted_graph_label.length();

      if (compacted_graph_label_length == 1) {
        // Add the neighbors of the chain to the neighbors of the last vertex.
        adjacency_list_[compacted_graph_vertex_id].insert(
            adjacency_list_[compacted_graph_vertex_id].end(),
            compacted_graph_adjacency_list_[compacted_graph_vertex_id].begin(),
            compacted_graph_adjacency_list_[compacted_graph_vertex_id].end());
      }

      for (GraphSizeType i = 1; i < compacted_graph_label_length; ++i) {
        labels_.emplace_back(compacted_graph_label[i]);
        adjacency_list_.emplace_back(std::vector<GraphSizeType>());

        // If this is the second vertex in the chain, add the link from the
        // first vertex to the second vertex in the chain.
        if (i == 1) {
          adjacency_list_[compacted_graph_vertex_id].push_back(vertex_id);
        }

        // If this is not the last vertex in the chain, which means it has a
        // next vertex, add link from current vertex to its next.
        if (i + 1 < compacted_graph_label_length) {
          adjacency_list_[vertex_id].push_back(vertex_id + 1);
        } else {
          // If this is the last vertex in the chain, add the neighbors of the
          // chain to the neighbors of the last vertex.
          adjacency_list_[vertex_id].insert(
              adjacency_list_[vertex_id].end(),
              compacted_graph_adjacency_list_[compacted_graph_vertex_id]
                  .begin(),
              compacted_graph_adjacency_list_[compacted_graph_vertex_id].end());
        }

        ++vertex_id;
      }

      ++compacted_graph_vertex_id;
    }
    // after the loop, compacted_graph_vertex_id should be the number of
    // vertices in the compacted graph
    // assert(compacted_graph_vertex_id == GetNumVerticesInCompactedGraph());
    // add an edge between the dummy and every other vertex
    // adjacency_list_[0].reserve(vertex_id - 1);
    // for (GraphSizeType i = 1; i < vertex_id; ++i) {
    //  adjacency_list_[0].push_back(i);
    //}
  }

  inline char GetComplementaryVertexLabelInCompactedGraph(
      GraphSizeType vertex_id, GraphSizeType vertex_j) const {
    return base_complement_[(int)compacted_graph_labels_[vertex_id][vertex_j]];
  }

  inline char GetReverseComplementaryVertexLabel(GraphSizeType vertex) const {
    return base_complement_[(int)labels_[vertex]];
  }

  inline char GetVertexLabel(GraphSizeType vertex) const {
    return labels_[vertex];
  }

  inline const std::vector<GraphSizeType> &GetNeighbors(
      GraphSizeType vertex) const {
    return adjacency_list_[vertex];
  }

  inline const std::vector<GraphSizeType> &GetNeighborsInCompatedGraph(
      GraphSizeType vertex) const {
    return compacted_graph_adjacency_list_[vertex];
  }

  void GenerateReverseComplementaryCharLabeledGraph() {
    reverse_complementary_adjacency_list_.assign(adjacency_list_.size(),
                                                 std::vector<GraphSizeType>());
    for (GraphSizeType vertex = 0; vertex < adjacency_list_.size(); ++vertex) {
      for (const GraphSizeType neighbor : adjacency_list_[vertex]) {
        reverse_complementary_adjacency_list_[neighbor].push_back(vertex);
      }
    }
  }

 protected:
  char base_complement_[256] = {
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 'T', 4, 'G', 4, 4, 4, 'C', 4,   4, 4, 4,
      4, 4, 'N', 4, 4,   4, 4, 4, 'A', 4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 'T', 4, 'G', 4, 4, 4, 'C', 4, 4,   4, 4, 4, 4,   'N', 4, 4, 4,
      4, 4, 'A', 4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4};
  int8_t base_to_int_[256] = {
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

  // For graph representation
  std::vector<std::vector<GraphSizeType>> adjacency_list_;
  std::vector<char> labels_;

  std::vector<std::vector<GraphSizeType>> compacted_graph_adjacency_list_;
  std::vector<std::string> compacted_graph_labels_;

  // For reverse complementary graph representation. We make the vertex ids
  // stable, and thus the vertex labels are just the corresponding reverse
  // complementary base in labels_. However, we construct a separate adjacency
  // list for fast neighbor queries in the reverse complementary graph.
  std::vector<std::vector<GraphSizeType>> reverse_complementary_adjacency_list_;
};

}  // namespace sgat
#endif  // SGAT_SEQUENCEGRAPH_H
