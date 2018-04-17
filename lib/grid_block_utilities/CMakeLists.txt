add_library(grid_block_utilities
  grid_block_utilities.f90
  test_block_utilities.f90
)
target_link_libraries(grid_block_utilities
    global
    grid_block_structure
    xyzq_io
)
target_include_directories(grid_block_utilities PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
