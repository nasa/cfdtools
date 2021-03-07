add_executable(gradient_based
  gradient_based.f90
  lcsfit2.f90
)
target_link_libraries(gradient_based cfdtools_global)
target_link_libraries(gradient_based
  cfdtools_numodules
  cfdtools_interplib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME gradient_based)
install(TARGETS gradient_based DESTINATION ${CMAKE_INSTALL_BINDIR})
