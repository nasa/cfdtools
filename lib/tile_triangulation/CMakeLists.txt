add_library(tile_triangulation
  triangulate_tiles.f90
)
target_link_libraries(tile_triangulation PRIVATE cfdtools_global)
target_include_directories(tile_triangulation PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
