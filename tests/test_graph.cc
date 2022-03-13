#include "dijkstra.h"
#include "gtest/gtest.h"
#include "navarro.h"
#include "scoring_schema.h"
#include "sequence_batch.h"
#include "sequence_graph.h"
#include "sga.h"

namespace sgat_testing {

class SequenceGraphTest : public ::testing::Test {
 protected:
  SequenceGraphTest()
      : txt_sequence_graph_(), gfa_sequence_graph_(), sequence_batch_(5) {}

  ~SequenceGraphTest() override {}

  void SetUp() override {
    sequence_batch_.InitializeLoading(sequence_file_path_);

    txt_sequence_graph_.LoadFromTxtFile(txt_sequence_graph_file_path_);
    txt_sequence_graph_.GenerateCharLabeledGraph();

    gfa_sequence_graph_.LoadFromGfaFile(gfa_sequence_graph_file_path_);
    gfa_sequence_graph_.GenerateCharLabeledGraph();
  }

  void TearDown() override { sequence_batch_.FinalizeLoading(); }

  sgat::SequenceGraph<int32_t> txt_sequence_graph_;
  const std::string txt_sequence_graph_file_path_ = "BRCA1_seq_graph.txt";

  sgat::SequenceGraph<int32_t> gfa_sequence_graph_;
  const std::string gfa_sequence_graph_file_path_ =
      "BRCA1_seq_graph_generated_from_txt.gfa";

  sgat::SequenceBatch sequence_batch_;
  const std::string sequence_file_path_ = "BRCA1_5_reads.fastq";
};

TEST_F(SequenceGraphTest, LoadFromTxtFileTest) {
  const uint32_t num_vertices =
      txt_sequence_graph_.GetNumVerticesInCompactedGraph() - 1;
  EXPECT_EQ(num_vertices, (uint32_t)2320)
      << "Number of vertices loaded is wrong! It should be 2320 but it is "
      << num_vertices;
  const uint32_t num_edges = txt_sequence_graph_.GetNumEdgesInCompactedGraph();
  EXPECT_EQ(num_edges, (uint32_t)2320)
      << "Number of edges loaded is wrong! It should be 2320 but it is "
      << num_edges;
}

TEST_F(SequenceGraphTest, LoadFromGfaFileTest) {
  const uint32_t num_vertices =
      gfa_sequence_graph_.GetNumVerticesInCompactedGraph() - 1;
  EXPECT_EQ(num_vertices, (uint32_t)2320)
      << "Number of vertices loaded is wrong! It should be 2320 but it is "
      << num_vertices;
  const uint32_t num_edges = gfa_sequence_graph_.GetNumEdgesInCompactedGraph();
  EXPECT_EQ(num_edges, (uint32_t)2320)
      << "Number of edges loaded is wrong! It should be 2320 but it is "
      << num_edges;
}

TEST_F(SequenceGraphTest, GenerateCharLabeledGraphOnTxtGraphTest) {
  const uint32_t num_vertices = txt_sequence_graph_.GetNumVertices() - 1;
  EXPECT_EQ(num_vertices, (uint32_t)139189)
      << "Number of vertices loaded is wrong! It should be 139189 but it is "
      << num_vertices;
  const uint32_t num_edges = txt_sequence_graph_.GetNumEdges();
  EXPECT_EQ(num_edges, (uint32_t)(139189))
      << "Number of edges loaded is wrong! It should be 139189 but it is "
      << num_edges;
}

TEST_F(SequenceGraphTest, GenerateCharLabeledGraphOnGfaGraphTest) {
  const uint32_t num_vertices = gfa_sequence_graph_.GetNumVertices() - 1;
  EXPECT_EQ(num_vertices, (uint32_t)139189)
      << "Number of vertices loaded is wrong! It should be 139189 but it is "
      << num_vertices;
  const uint32_t num_edges = gfa_sequence_graph_.GetNumEdges();
  EXPECT_EQ(num_edges, (uint32_t)(139189))
      << "Number of edges loaded is wrong! It should be 139189 but it is "
      << num_edges;
}

TEST_F(SequenceGraphTest, AlignUsingLinearGapPenaltyTest) {
  const uint32_t num_loaded_sequences = sequence_batch_.LoadBatch();
  const int32_t max_alignment_scores[5] = {62, 25, 54, 9, 37};

  const sgat::ScoringSchema<int32_t> scoring_schema;

  sgat::SgaAligner</*GraphSizeType=*/int32_t, /*QueryLengthType=*/int32_t,
                   /*ScoreType=*/int32_t>
      sga_aligner;

  for (uint32_t i = 0; i < num_loaded_sequences; ++i) {
    const int32_t alignment_score = sga_aligner.AlignUsingLinearGapPenalty(
        i, sequence_batch_, txt_sequence_graph_, scoring_schema);
    EXPECT_EQ(alignment_score, max_alignment_scores[i])
        << "Alignment score for sequence" << i << " is wrong! It should be "
        << max_alignment_scores[i] << " but it is " << alignment_score;
  }
}

TEST_F(SequenceGraphTest, AlignUsingLinearGapPenaltyWithNavarroAlgorithmTest) {
  const uint32_t num_loaded_sequences = sequence_batch_.LoadBatch();
  const int32_t max_alignment_scores[5] = {62, 25, 54, 9, 37};
  const sgat::ScoringSchema<int32_t> scoring_schema;

  sgat::NavarroAligner</*GraphSizeType=*/int32_t, /*QueryLengthType=*/int32_t,
                       /*ScoreType=*/int32_t>
      navarro_aligner;

  for (uint32_t i = 0; i < num_loaded_sequences; ++i) {
    const int32_t alignment_score =
        navarro_aligner
            .AlignUsingLinearGapPenaltyWithNavarroAlgorithmSemiGloballyOnTwoDirections(
                i, sequence_batch_, txt_sequence_graph_, scoring_schema);
    EXPECT_EQ(alignment_score, max_alignment_scores[i])
        << "Alignment score for sequence" << i << " is wrong! It should be "
        << max_alignment_scores[i] << " but it is " << alignment_score;
  }
}

TEST_F(SequenceGraphTest,
       AlignUsingLinearGapPenaltyWithNavarroAlgorithmOnGfaGraphTest) {
  const uint32_t num_loaded_sequences = sequence_batch_.LoadBatch();
  const int32_t max_alignment_scores[5] = {62, 25, 54, 9, 37};
  const sgat::ScoringSchema<int32_t> scoring_schema;

  sgat::NavarroAligner</*GraphSizeType=*/int32_t, /*QueryLengthType=*/int32_t,
                       /*ScoreType=*/int32_t>
      navarro_aligner;

  for (uint32_t i = 0; i < num_loaded_sequences; ++i) {
    const int32_t alignment_score =
        navarro_aligner
            .AlignUsingLinearGapPenaltyWithNavarroAlgorithmSemiGloballyOnTwoDirections(
                i, sequence_batch_, gfa_sequence_graph_, scoring_schema);
    EXPECT_EQ(alignment_score, max_alignment_scores[i])
        << "Alignment score for sequence" << i << " is wrong! It should be "
        << max_alignment_scores[i] << " but it is " << alignment_score;
  }
}

TEST_F(SequenceGraphTest, AlignUsingLinearGapPenaltyWithDijkstraAlgorithmTest) {
  const uint32_t num_loaded_sequences = sequence_batch_.LoadBatch();
  const int32_t max_alignment_scores[5] = {62, 25, 54, 9, 37};

  const sgat::ScoringSchema<int32_t> scoring_schema;

  sgat::DijkstraAligner</*GraphSizeType=*/int32_t, /*QueryLengthType=*/int32_t,
                        /*ScoreType=*/int32_t>
      dijkstra_aligner;

  for (uint32_t i = 0; i < num_loaded_sequences; ++i) {
    const int32_t alignment_score =
        dijkstra_aligner.AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
            i, sequence_batch_, txt_sequence_graph_, scoring_schema);
    EXPECT_EQ(alignment_score, max_alignment_scores[i])
        << "Alignment score for sequence" << i << " is wrong! It should be "
        << max_alignment_scores[i] << " but it is " << alignment_score;
  }
}
}  // namespace sgat_testing

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
