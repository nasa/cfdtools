add_executable(combine
  combine.f
)
target_link_libraries(combine
  global
  interplib
  searchlib
  numodules
)
install(TARGETS combine DESTINATION ${CMAKE_INSTALL_BINDIR})
