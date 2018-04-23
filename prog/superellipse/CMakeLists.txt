add_executable(superellipse
  superellipse.f
)
target_link_libraries(superellipse global)
install(TARGETS superellipse DESTINATION ${CMAKE_INSTALL_BINDIR})
