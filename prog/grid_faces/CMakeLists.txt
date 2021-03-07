cfdtools_add_executable(grid_faces
  grid_faces.f90
)
target_link_libraries(cfdtools_grid_faces
  cfdtools_xyzq_io
)
