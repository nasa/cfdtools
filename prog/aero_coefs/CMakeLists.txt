add_executable(aero_coefs
    aero_coefs.f90
    patch_coefs.f90
)

target_link_libraries(aero_coefs
  PRIVATE
    cfdtools_global
    cfdtools_gridlib
    cfdtools_geomlib
    cfdtools_progtools
    cfdtools_xyzq_io
)

install(TARGETS aero_coefs DESTINATION ${CMAKE_INSTALL_BINDIR})

