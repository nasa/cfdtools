add_executable(surface_uvs
  surface_uvs.f90
)
target_link_libraries(surface_uvs cfdtools_global cfdtools_geomlib cfdtools_xyzq_io)
install(TARGETS surface_uvs DESTINATION ${CMAKE_INSTALL_BINDIR})
