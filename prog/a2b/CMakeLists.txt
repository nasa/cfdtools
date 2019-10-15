add_executable(a2b a2b.f90)
target_link_libraries(a2b PRIVATE global progtools)
install(TARGETS a2b DESTINATION ${CMAKE_INSTALL_BINDIR})

