add_library(ugridlib
  calcphi.f.unused
  cut_triangles.f
  cutorder.f
  cutorder_nf.f
  cutsurf.f
  cutsurf_nf.f
  cutsurfi.f
  cutsurfi_nf.f
  isintet.f.unused
  nonlinear_tri_interp.f90
  points_delaunay_naive_2d.f90
  search_tri_2d.f90
  solid_angle_quadrature_3.f90
  solid_angle_tri.f90
  solid_angle_tri.f90.single_rbar
  spherical_triangulation.f90
  testtet.f
  tet_coefs.f
  tet_volumes.f90
  tri_areas.f90
  tri_centers_to_vertices.f90
  tri_coefs.f
  triangulate_patches.f
)
target_link_libraries(ugridlib global)
install(FILES README DESTINATION doc)
