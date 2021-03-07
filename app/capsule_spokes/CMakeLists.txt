cfdtools_add_executable(capsule_spokes capsule_spokes.f90)
target_link_libraries(cfdtools_capsule_spokes
    cfdtools_geomlib
    cfdtools_numodules
    cfdtools_progtools
    cfdtools_table_io
)
