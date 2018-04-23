add_executable(compare_flows
  compare_flows.f90
)
target_link_libraries(compare_flows global)
install(FILES README DESTINATION doc RENAME compare_flows)
install(TARGETS compare_flows DESTINATION ${CMAKE_INSTALL_BINDIR})
