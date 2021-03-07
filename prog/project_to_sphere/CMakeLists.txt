add_executable(project_to_sphere
  project_to_sphere.f90
)
target_link_libraries(project_to_sphere cfdtools_global cfdtools_xyzq_io)
install(TARGETS project_to_sphere DESTINATION ${CMAKE_INSTALL_BINDIR})
