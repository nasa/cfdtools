add_executable(inertial_and_relative
  inertial_and_relative.f
)
target_link_libraries(inertial_and_relative global)
install(TARGETS inertial_and_relative DESTINATION ${CMAKE_INSTALL_BINDIR})
