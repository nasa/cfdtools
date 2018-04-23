add_executable(fourdigit
  fourdigit.f
)
target_link_libraries(fourdigit global)
install(TARGETS fourdigit DESTINATION ${CMAKE_INSTALL_BINDIR})
