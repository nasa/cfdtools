add_executable(npopt_driver
  confun.f90
  fun.f90
  npopt_driver.f90
  objfun.f90
)
target_link_libraries(npopt_driver cfdtools_global)
install(FILES README DESTINATION doc RENAME npopt_driver)
install(TARGETS npopt_driver DESTINATION ${CMAKE_INSTALL_BINDIR})
