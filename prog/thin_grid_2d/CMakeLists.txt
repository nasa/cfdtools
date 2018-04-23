add_executable(thin_grid_2d
  thin_grid_2d.f90
)
target_link_libraries(thin_grid_2d global)
target_link_libraries(thin_grid_2d
  xyzq_io
)
install(FILES README DESTINATION doc RENAME thin_grid_2d)
install(TARGETS thin_grid_2d DESTINATION ${CMAKE_INSTALL_BINDIR})
