add_executable(gu
  gu.f90
)
target_link_libraries(gu cfdtools_global)
install(FILES README DESTINATION doc RENAME gu)
install(TARGETS gu DESTINATION ${CMAKE_INSTALL_BINDIR})
