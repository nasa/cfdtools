cfdtools_add_executable(wedge wedge.f90)
target_link_libraries(cfdtools_wedge
    cfdtools_geomlib
    cfdtools_xyzq_io
)
