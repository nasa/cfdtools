cfdtools_add_executable(aero_coefs
    aero_coefs.f90
    patch_coefs.f90
)
target_link_libraries(cfdtools_aero_coefs
    cfdtools_gridlib
    cfdtools_geomlib
    cfdtools_progtools
    cfdtools_xyzq_io
)

