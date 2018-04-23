add_executable(surface_project
  surface_project.f90
)
target_link_libraries(surface_project global)
install(TARGETS surface_project DESTINATION ${CMAKE_INSTALL_BINDIR})
