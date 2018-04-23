add_executable(prepare_local_analysis
  reader.f
  upcase.f
  prepare_local_analysis.f90
)
target_link_libraries(prepare_local_analysis global)
target_link_libraries(prepare_local_analysis
  progtools
)
install(FILES README DESTINATION doc RENAME prepare_local_analysis)
install(TARGETS prepare_local_analysis DESTINATION ${CMAKE_INSTALL_BINDIR})
