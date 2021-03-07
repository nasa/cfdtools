cfdtools_add_library(triangulation_io
  triangulation_io.f90
  # TODO: Refactor tri_areas, tri_center_to_vertices into a new lib
  #       As written, there is a circular dependency between triangulation_io
  #       (which need the routines below) and ugridlib/spherical_triangulation.f90
  #       which requries the triangulation_io library. As a hack, we simply
  #       re-complile the utilities needed to triangulation_io to break the cycle
  ../ugridlib/tri_areas.f90
  ../ugridlib/tri_centers_to_vertices.f90
)
target_link_libraries(cfdtools_triangulation_io
  PRIVATE
    cfdtools_eigenlib
    cfdtools_progtools
)
