add_executable(convertq
  convertq.f
)
target_link_libraries(convertq global interplib searchlib)
install(TARGETS convertq DESTINATION ${CMAKE_INSTALL_BINDIR})
