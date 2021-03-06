add_executable(deconstruct
  deconstruct.f90
)
target_link_libraries(deconstruct cfdtools_global)
target_link_libraries(deconstruct
  xyzq_io
)
install(FILES README DESTINATION doc RENAME deconstruct)
install(TARGETS deconstruct DESTINATION ${CMAKE_INSTALL_BINDIR})
