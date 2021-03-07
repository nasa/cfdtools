add_executable(surface_peaks
  surface_peaks.f90
)
target_link_libraries(surface_peaks cfdtools_global cfdtools_tecplot_io cfdtools_progtools)
install(TARGETS surface_peaks DESTINATION ${CMAKE_INSTALL_BINDIR})
