add_executable(neqair_data
  neqair_data.f90
)
target_link_libraries(neqair_data cfdtools_global)
target_link_libraries(neqair_data
  cfdtools_xyzq_io
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_numodules
  cfdtools_progtools
)
install(TARGETS neqair_data DESTINATION ${CMAKE_INSTALL_BINDIR})
