add_executable(prepare_fiat_data
  prepare_fiat_data.f90
)
target_link_libraries(prepare_fiat_data global)
target_link_libraries(prepare_fiat_data
  table_io
  progtools
)
install(TARGETS prepare_fiat_data DESTINATION ${CMAKE_INSTALL_BINDIR})
