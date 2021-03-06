add_executable(revolve_grid
  revolve_grid.f90
)
target_link_libraries(revolve_grid cfdtools_global)
target_link_libraries(revolve_grid
  geomlib
  gridlib
  xyzq_io
)
install(FILES README DESTINATION doc RENAME revolve_grid)
install(TARGETS revolve_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
