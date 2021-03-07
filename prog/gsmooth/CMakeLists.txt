cfdtools_add_executable(gsmooth
  gsmooth.f90
)
target_link_libraries(cfdtools_gsmooth
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_xyzq_io
)
