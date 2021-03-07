add_executable(gsmooth
  gsmooth.f90
)
target_link_libraries(gsmooth cfdtools_global)
target_link_libraries(gsmooth
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_xyzq_io
)
install(FILES README DESTINATION doc RENAME gsmooth)
install(TARGETS gsmooth DESTINATION ${CMAKE_INSTALL_BINDIR})
