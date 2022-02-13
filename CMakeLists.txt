cmake_minimum_required(VERSION 3.3)

project(sgat DESCRIPTION "Sequence to graph alignment tools")

# sgat and GoogleTest require at least C++11
set(CMAKE_CXX_STANDARD 11)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

add_subdirectory(src)
add_subdirectory(apps)

option(SGAT_BUILD_TESTING "Build the tests")
if(SGAT_BUILD_TESTING)
  message(STATUS "Build the tests.")

  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  add_subdirectory(tests)
endif()

option(SGAT_BUILD_BENCHMARK "Build the benchmark")
if(SGAT_BUILD_BENCHMARK)
  message(STATUS "Build the benchmark.")

  add_subdirectory(extern/benchmark)
  add_subdirectory(benchmark)
endif()
