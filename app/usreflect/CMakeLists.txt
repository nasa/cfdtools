cfdtools_add_executable(uslos usreflect.f90)
target_link_libraries(cfdtools_reflect
    cfdtools_triangulation_io
    cfdtools_ugridlib
    cfdtools_geomlib
    cfdtools_eigenlib
    cfdtools_progtools
)
