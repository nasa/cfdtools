add_executable(prepare_neqair_data
  prepare_neqair_data.f90
)
target_link_libraries(prepare_neqair_data cfdtools_global)
target_link_libraries(prepare_neqair_data
  gridlib
  progtools
)
install(FILES README DESTINATION doc RENAME prepare_neqair_data)
install(TARGETS prepare_neqair_data DESTINATION ${CMAKE_INSTALL_BINDIR})
