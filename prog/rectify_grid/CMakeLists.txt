add_executable(rectify_grid
  rectify_grid.f90
)
target_link_libraries(rectify_grid cfdtools_global)
target_link_libraries(rectify_grid
  gridlib
  progtools
  geomlib
  xyzq_io
)
install(FILES README DESTINATION doc RENAME rectify_grid)
install(TARGETS rectify_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
