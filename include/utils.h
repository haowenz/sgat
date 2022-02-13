#ifndef SGAT_UTILS_H_
#define SGAT_UTILS_H_

#include <sys/resource.h>
#include <sys/time.h>

#include <cstdint>
#include <iostream>
#include <vector>

namespace sgat {

// For timing
inline static double GetRealTime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return tp.tv_sec + tp.tv_usec * 1e-6;
}

inline static double GetCPUTime() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
         1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

inline static void ExitWithMessage(const std::string &message) {
  std::cerr << message << std::endl;
  exit(-1);
}

}  // namespace sgat

#endif  // SGAT_UTILS_H_
