add_library(ugridlib
  cut_triangles.f
  cutorder.f
  cutorder_nf.f
  cutsurf.f
  cutsurf_nf.f
  cutsurfi.f
  cutsurfi_nf.f
  nonlinear_tri_interp.f90
  points_delaunay_naive_2d.f90
  search_tri_2d.f90
  solid_angle_quadrature_3.f90
  solid_angle_tri.f90
  spherical_triangulation.f90
  testtet.f
  tet_coefs.f
  tet_volumes.f90
  tri_areas.f90
  tri_centers_to_vertices.f90
  tri_coefs.f
  triangulate_patches.f
)
target_link_libraries(ugridlib global triangulation_io)
target_include_directories(ugridlib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc)
