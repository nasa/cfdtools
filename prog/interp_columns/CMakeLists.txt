add_executable(interp_columns
  interp_columns.f90
)
target_link_libraries(interp_columns cfdtools_global)
target_link_libraries(interp_columns
  cfdtools_numodules
  cfdtools_interplib
  cfdtools_progtools
  cfdtools_linsys
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME interp_columns)
install(TARGETS interp_columns DESTINATION ${CMAKE_INSTALL_BINDIR})
