#include <string>

#include "dijkstra.h"
#include "sequence_batch.h"
#include "sequence_graph.h"
#include "utils.h"

int main(int argc, char *argv[]) {
  std::string sequence_graph_file_path;
  std::string sequence_file_path;
  if (argc != 3) {
    std::cerr << "Usage:\t" << argv[0] << "\tgraph_file\tread_file\n";
    exit(-1);
  } else {
    sequence_graph_file_path = argv[1];
    sequence_file_path = argv[2];
  }

  sgat::SequenceGraph<int32_t> sequence_graph;
  sequence_graph.LoadFromGfaFile(sequence_graph_file_path);
  sequence_graph.GenerateCharLabeledGraph();

  const uint32_t max_batch_size = 1000000;
  sgat::SequenceBatch sequence_batch(max_batch_size);
  sequence_batch.InitializeLoading(sequence_file_path);

  const sgat::ScoringSchema<int32_t> scoring_schema;
  sgat::DijkstraAligner<int32_t, int32_t, int32_t> dijkstra_aligner;

  uint32_t num_sequences = sequence_batch.LoadBatch();
  uint64_t num_total_sequences = 0;

  double mapping_start_real_time = sgat::GetRealTime();

  while (num_sequences > 0) {
    for (uint32_t si = 0; si < num_sequences; ++si) {
      dijkstra_aligner.AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
          si, sequence_batch, sequence_graph, scoring_schema);
    }
    num_total_sequences += num_sequences;
    num_sequences = sequence_batch.LoadBatch();
  }

  std::cerr << "Mapped " << num_total_sequences << " sequences in "
            << sgat::GetRealTime() - mapping_start_real_time << "s"
            << std::endl;

  sequence_batch.FinalizeLoading();
}
