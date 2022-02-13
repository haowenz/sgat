find_package(ZLIB REQUIRED)

include(ExternalProject)
ExternalProject_Add(gfatools
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/extern/gfatools"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${CMAKE_COMMAND} -E env
  HOST_CC=${CMAKE_CXX_COMPILER}
  CC=${CMAKE_CXX_COMPILER}
  make -C <SOURCE_DIR> libgfa1.a
  BUILD_IN_SOURCE TRUE
  INSTALL_COMMAND ""
  )

add_library(sgat STATIC sequence_batch.cc)
target_include_directories(sgat PUBLIC "${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/extern/klib" "${CMAKE_SOURCE_DIR}/extern/gfatools")
target_link_libraries(sgat ZLIB::ZLIB "${CMAKE_SOURCE_DIR}/extern/gfatools/libgfa1.a")
add_dependencies(sgat gfatools)
