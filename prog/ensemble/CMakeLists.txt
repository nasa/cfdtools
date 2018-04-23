add_executable(ensemble
  ensemble.f90
  rotation_matrix.f90
)
target_link_libraries(ensemble global)
install(TARGETS ensemble DESTINATION ${CMAKE_INSTALL_BINDIR})
