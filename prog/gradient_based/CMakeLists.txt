add_executable(gradient_based
  gradient_based.f90
  lcsfit2.f90
)
target_link_libraries(gradient_based global)
target_link_libraries(gradient_based
  numodules
  interplib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME gradient_based)
install(TARGETS gradient_based DESTINATION ${CMAKE_INSTALL_BINDIR})
