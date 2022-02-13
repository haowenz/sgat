#include <string>

#include "sequence_graph.h"

int main(int argc, char *argv[]) {
  std::string sequence_graph_file_path;
  std::string gfa_sequence_graph_output_file_path;
  if (argc != 3) {
    std::cerr << "Usage:\t" << argv[0] << "\tgraph_file\toutput_file\n";
    exit(-1);
  } else {
    sequence_graph_file_path = argv[1];
    gfa_sequence_graph_output_file_path = argv[2];
  }
  sgat::SequenceGraph<int32_t> sequence_graph;
  sequence_graph.LoadFromGfaFile(sequence_graph_file_path);
  sequence_graph.GenerateCharLabeledGraph();

  sequence_graph.OutputCharLabeledGraphInGFA(
      gfa_sequence_graph_output_file_path);
}
