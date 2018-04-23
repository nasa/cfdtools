add_executable(peaker
  peaker.f
)
target_link_libraries(peaker global)
install(TARGETS peaker DESTINATION ${CMAKE_INSTALL_BINDIR})
