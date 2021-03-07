add_executable(grid_faces
  grid_faces.f90
)
target_link_libraries(grid_faces cfdtools_global)
target_link_libraries(grid_faces
  cfdtools_xyzq_io
)
install(FILES README DESTINATION doc RENAME grid_faces)
install(TARGETS grid_faces DESTINATION ${CMAKE_INSTALL_BINDIR})
