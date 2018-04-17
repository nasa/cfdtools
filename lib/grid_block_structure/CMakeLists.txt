add_library(grid_block_structure
  grid_block_structure.f90
)
target_link_libraries(grid_block_structure global)
target_include_directories(grid_block_structure PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
