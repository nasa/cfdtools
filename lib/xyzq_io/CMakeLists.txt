add_library(xyzq_io
  f_io.f90
  xyq_io.f90
  xyzq_io.f90
)
target_link_libraries(xyzq_io
  PRIVATE
    global
  PUBLIC
    grid_block_structure
)
target_include_directories(xyzq_io PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME xyzq_io)
