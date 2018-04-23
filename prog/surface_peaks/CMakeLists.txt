add_executable(surface_peaks
  surface_peaks.f90
)
target_link_libraries(surface_peaks global)
install(TARGETS surface_peaks DESTINATION ${CMAKE_INSTALL_BINDIR})
