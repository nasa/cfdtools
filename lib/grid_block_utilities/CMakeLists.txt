add_library(grid_block_utilities
  grid_block_utilities.f90
)
target_link_libraries(grid_block_utilities
  PRIVATE
    global
    xyzq_io
  PUBLIC
    grid_block_structure
)
target_include_directories(grid_block_utilities PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
