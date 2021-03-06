add_executable(nose_radius
  convert_variables.f
  nose_radius.f
)
target_link_libraries(nose_radius cfdtools_global)
install(TARGETS nose_radius DESTINATION ${CMAKE_INSTALL_BINDIR})
