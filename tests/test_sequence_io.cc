#include "gtest/gtest.h"
#include "sequence_batch.h"

namespace sgat_testing {

class SequenceBatchTest : public ::testing::Test {
 protected:
  SequenceBatchTest() : sequence_batch_(2) {
    sequence_batch_.InitializeLoading(sequence_file_path_);
  }

  ~SequenceBatchTest() override { sequence_batch_.FinalizeLoading(); }

  void SetUp() override {}

  void TearDown() override {}

  sgat::SequenceBatch sequence_batch_;
  const std::string sequence_file_path_ = "BRCA1_5_reads.fastq";
};

TEST_F(SequenceBatchTest, LoadBatchTest) {
  uint32_t num_loaded_sequences = sequence_batch_.LoadBatch();
  ASSERT_EQ(num_loaded_sequences, (uint32_t)2)
      << "Number of sequences loaded is wrong! It should be 2 but it is "
      << num_loaded_sequences;
  num_loaded_sequences = sequence_batch_.LoadBatch();
  ASSERT_EQ(num_loaded_sequences, (uint32_t)2)
      << "Number of sequences loaded is wrong! It should be 2 but it is "
      << num_loaded_sequences;
  num_loaded_sequences = sequence_batch_.LoadBatch();
  ASSERT_EQ(num_loaded_sequences, (uint32_t)1)
      << "Number of sequences loaded is wrong! It should be 2 but it is "
      << num_loaded_sequences;

  const char *true_sequence =
      "TCTTGTAATTTAATTTCGATTACTAATTTCCTGAAAATTATAGAACTAGATAAAGCTATATAGTTGTGGATT"
      "ATTTATGGTATAGTTTACTTGAGAAATAATTATTAAATATTAGTGGAAAAAGCTATACTTTGGGTATGATAA"
      "GGAACTTCCTCAATTGAATTTCCTTTCCTATCTGTAAAAAGCAAGTAGGGTTAAATAGTTTTATTCCGCCAG"
      "AAGGCATCTTTTATTCTTCTCCCCTTGTCCTCACATGGGTGAATTTACCAGCACATTTAACTAAATTCAGAA"
      "ACTGGTTCCAAATGTACTGCAGATAGTAGGCAATTTC";

  const char *true_sequence_qual =
      ")**)**)********)*******)**********)*)****************)*******)***)******"
      "*****************)*)***)*)***)*****)****************************)*******"
      "************)****)*************)******)*))*)*******)*)***))***)****)****"
      "*)**)*******))*********)**))******)***))****)**************************)"
      "************************)*****)****)*";

  const char *loaded_sequence = sequence_batch_.GetSequenceAt(0);
  const char *loaded_sequence_qual = sequence_batch_.GetSequenceQualAt(0);

  ASSERT_STREQ(loaded_sequence, true_sequence);
  ASSERT_STREQ(loaded_sequence_qual, true_sequence_qual);
}

}  // namespace sgat_testing

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
