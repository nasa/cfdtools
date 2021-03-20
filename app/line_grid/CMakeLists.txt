cfdtools_add_executable(line_grid
  line_grid.f90
)
target_link_libraries(cfdtools_line_grid cfdtools_gridlib)
