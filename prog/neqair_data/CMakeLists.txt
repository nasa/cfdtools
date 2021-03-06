add_executable(neqair_data
  neqair_data.f90
)
target_link_libraries(neqair_data cfdtools_global)
target_link_libraries(neqair_data
  xyzq_io
  gridlib
  geomlib
  interplib
  searchlib
  numodules
  progtools
)
install(TARGETS neqair_data DESTINATION ${CMAKE_INSTALL_BINDIR})
