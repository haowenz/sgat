add_executable(SequenceIOTest test_sequence_io.cc)
target_link_libraries(SequenceIOTest sgat gtest_main)

add_executable(SequenceGraphTest test_graph.cc)
target_link_libraries(SequenceGraphTest sgat gtest_main)

include(GoogleTest)
message(STATUS "cmake src dir: ${CMAKE_SOURCE_DIR}")
gtest_discover_tests(SequenceIOTest WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/data")
gtest_discover_tests(SequenceGraphTest WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/data")
