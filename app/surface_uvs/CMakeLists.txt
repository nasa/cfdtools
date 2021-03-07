cfdtools_add_executable(surface_uvs
  surface_uvs.f90
)
target_link_libraries(cfdtools_surface_uvs
  cfdtools_geomlib
  cfdtools_xyzq_io
)
