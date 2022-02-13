#include <string>

#include "navarro.h"
#include "scoring_schema.h"
#include "sequence_batch.h"
#include "sequence_graph.h"
#include "utils.h"

int main(int argc, char *argv[]) {
  std::string sequence_graph_file_path;
  std::string sequence_file_path;
  int32_t start_vertex = 1;
  if (argc != 4) {
    std::cerr << "Usage:\t" << argv[0]
              << "\t1-based_start_vertex_id\tgraph_file\tread_file\n";
    exit(-1);
  } else {
    start_vertex = atoi(argv[1]);
    sequence_graph_file_path = argv[2];
    sequence_file_path = argv[3];
  }
  sgat::SequenceGraph<int32_t> sequence_graph;
  sequence_graph.LoadFromGfaFile(sequence_graph_file_path);
  sequence_graph.GenerateCharLabeledGraph();

  const uint32_t max_batch_size = 1000000;
  sgat::SequenceBatch sequence_batch(max_batch_size);
  sequence_batch.InitializeLoading(sequence_file_path);

  const sgat::ScoringSchema<int32_t> scoring_schema;
  sgat::NavarroAligner<int32_t, int32_t, int32_t> navarro_aligner;

  uint32_t num_sequences = sequence_batch.LoadBatch();
  uint64_t num_total_sequences = 0;

  double mapping_start_real_time = sgat::GetRealTime();

  while (num_sequences > 0) {
    for (uint32_t si = 0; si < num_sequences; ++si) {
      navarro_aligner.ForwardExtendUsingLinearGapPenaltyWithNavarroAlgorithm(
          si, sequence_batch, sequence_graph, start_vertex, scoring_schema);
    }
    num_total_sequences += num_sequences;
    num_sequences = sequence_batch.LoadBatch();
  }

  std::cerr << "Extend " << num_total_sequences << " sequences in "
            << sgat::GetRealTime() - mapping_start_real_time
            << "s, starting from vertex " << start_vertex << std::endl;

  sequence_batch.FinalizeLoading();
}
