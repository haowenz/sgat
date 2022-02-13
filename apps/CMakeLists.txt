add_executable(sgat_example example.cc)
target_link_libraries(sgat_example PRIVATE sgat)

add_executable(txt2gfa txt2gfa.cc)
target_link_libraries(txt2gfa PRIVATE sgat)

add_executable(gfa2char_gfa gfa2char_gfa.cc)
target_link_libraries(gfa2char_gfa PRIVATE sgat)

add_executable(dijkstra_extend dijkstra_extend.cc)
target_link_libraries(dijkstra_extend PRIVATE sgat)

add_executable(navarro_extend navarro_extend.cc)
target_link_libraries(navarro_extend PRIVATE sgat)
