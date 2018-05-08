# This library is a copy of the gridlib library but compiled to support the
# extended versions of grid_block_structure defined in the tecplot_io library.
# In a perfect world, the tecplot version of these types would be extensions of
# the base type and the gridlib library would work with either version.
# However, the ability to do this sort generic programming is a recent addition
# to Fortran, so we would have to do an extensive rewrite of gridlib and client
# libraries to do it. Instead, we just compile mulitple versions of the
# library, since this is how CFDTOOLS has historicall solved the problem.

add_library(tecplot_gridlib
  # Use same sources grid_lib
  ../gridlib/adjustn.f90
  ../gridlib/adjustn2.f90
  ../gridlib/arbdis.f
  ../gridlib/arcdis.f
  ../gridlib/average_increments.f90
  ../gridlib/blgrid.f
  ../gridlib/catenary_grid.f90
  ../gridlib/cellvol.f
  ../gridlib/changen.f
  ../gridlib/changen1d.f
  ../gridlib/changen2d.f
  ../gridlib/compass.f
  ../gridlib/condis.f
  ../gridlib/copy_face.f90
  ../gridlib/curvdis.f
  ../gridlib/curvdis2.f90
  ../gridlib/curvdis3d.f
  ../gridlib/curvdis3d2.f90
  ../gridlib/curvdis_behavior.f90
  ../gridlib/curvdis_l.f
  ../gridlib/delq3d.f
  ../gridlib/determine_grid_dim.f90
  ../gridlib/determine_grid_form.f90
  ../gridlib/distrib.f
  ../gridlib/dstrib.f
  ../gridlib/ecompass.f
  ../gridlib/ellip2d.f
  ../gridlib/ellip3d.f
  ../gridlib/ellipq3d.f
  ../gridlib/expdis.f
  ../gridlib/expdis2.f
  ../gridlib/expdis4.f
  ../gridlib/expdis5.f
  ../gridlib/fixogrid.f
  ../gridlib/flow_gradients.f90
  ../gridlib/flow_gradients_nf.f90
  ../gridlib/foilgrd.f
  ../gridlib/foilgrid.f
  ../gridlib/geodis.f
  ../gridlib/geodis2.f
  ../gridlib/graddis.f
  ../gridlib/graddis2.f90
  ../gridlib/graddis3d.f
  ../gridlib/graddis3d2.f90
  ../gridlib/gridspec.f
  ../gridlib/htdis2.f
  ../gridlib/htdis4.f
  ../gridlib/lapgrid.f
  ../gridlib/min_max_curv_edge_2d.f90
  ../gridlib/modgrid2.f
  ../gridlib/modgrid3.f
  ../gridlib/modgrid4.f
  ../gridlib/morph_line_3d.f
  ../gridlib/nuline2d.f
  ../gridlib/nuline3d.f
  ../gridlib/permute_block.f90
  ../gridlib/permute_block_2d.f90
  ../gridlib/shockgrid.f
  ../gridlib/skew2d.f
  ../gridlib/smooth1d.f
  ../gridlib/smoothx.f
  ../gridlib/stretch.f90
  ../gridlib/tfi2d.f
  ../gridlib/tfi3d.f
  ../gridlib/tfint3d.f
  ../gridlib/tfint3f.f
  ../gridlib/tfiq3d.f
  ../gridlib/tfiq3xyz.f
  ../gridlib/uniform_edge.f90
  ../gridlib/vertices_to_centers.f90
  ../gridlib/vertices_to_centers_2d.f90
  ../gridlib/vinokur.f
  ../gridlib/warp2d.f
  ../gridlib/warp3d.f
  ../gridlib/warpq3d.f
  ../gridlib/warpq3d2.f90
  ../gridlib/xgrid.f
)
target_link_libraries(tecplot_gridlib
  PRIVATE
    global
    intrinsics
  PUBLIC
    # Link against Tecplot-customized verions of grid_block_structure
    tecplot_io
)
target_include_directories(tecplot_gridlib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
