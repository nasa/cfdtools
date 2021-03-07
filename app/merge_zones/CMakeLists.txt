cfdtools_add_executable(merge_zones merge_zones.f90)
target_link_libraries(cfdtools_merge_zones
    cfdtools_progtools
    cfdtools_tecplot_io
)
