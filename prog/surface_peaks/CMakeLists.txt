cfdtools_add_executable(surface_peaks
  surface_peaks.f90
)
target_link_libraries(cfdtools_surface_peaks
  cfdtools_tecplot_io
  cfdtools_progtools
)
