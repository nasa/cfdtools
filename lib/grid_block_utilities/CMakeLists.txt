add_library(grid_block_utilities
  grid_block_structure.f90
  grid_block_utilities.f90
  test_block_utilities.f90
)
target_link_libraries(grid_block_utilities global)
install(FILES README DESTINATION doc)
