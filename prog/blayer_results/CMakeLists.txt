add_executable(blayer_results
  blayer_results.f90
)
target_link_libraries(blayer_results global)
install(TARGETS blayer_results DESTINATION ${CMAKE_INSTALL_BINDIR})
