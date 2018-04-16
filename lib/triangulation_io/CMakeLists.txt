add_library(triangulation_io
  triangulation_io.f90
  triangulation_io.f90.before-area-volume-utilities
  triangulation_io.f90.before-volume-io
  triangulation_io.f90.before-ZONETYPE=FETriangle
  triangulation_io.f90.before.NASTRAN.options
)
target_link_libraries(triangulation_io global)
install(FILES README DESTINATION doc)
