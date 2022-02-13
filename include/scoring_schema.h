#ifndef SGAT_SCORINGSCHEMA_H
#define SGAT_SCORINGSCHEMA_H

#include <cstdint>

namespace sgat {

template <class ScoreType = int16_t>
struct ScoringSchema {
  ScoreType substitution_penalty = 1;
  ScoreType deletion_penalty = 1;
  ScoreType insertion_penalty = 1;
};

}  // namespace sgat
#endif  // SGAT_SCORINGSCHEMA_H
