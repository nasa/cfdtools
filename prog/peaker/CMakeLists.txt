add_executable(peaker
  peaker.f
)
target_link_libraries(peaker cfdtools_global)
install(TARGETS peaker DESTINATION ${CMAKE_INSTALL_BINDIR})
