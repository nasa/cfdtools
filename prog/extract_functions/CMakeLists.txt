add_executable(extract_functions
  extract_functions.f90
)
target_link_libraries(extract_functions cfdtools_global)
target_link_libraries(extract_functions
  cfdtools_progtools
  cfdtools_gridlib
  cfdtools_xyzq_io
)
install(TARGETS extract_functions DESTINATION ${CMAKE_INSTALL_BINDIR})
