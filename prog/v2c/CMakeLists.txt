add_executable(v2c
  v2c.f90
)
target_link_libraries(v2c cfdtools_global)
target_link_libraries(v2c
  progtools
  gridlib
  xyzq_io
)
install(FILES README DESTINATION doc RENAME v2c)
install(TARGETS v2c DESTINATION ${CMAKE_INSTALL_BINDIR})
