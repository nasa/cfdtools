add_executable(ensemble
  ensemble.f90
)
target_link_libraries(ensemble cfdtools_global)
target_link_libraries(ensemble
    intrinsics
    tecplot_adt
    tecplot_gridlib
    tecplot_io
    tecplot_surface_patch_utilities
    progtools
)
install(TARGETS ensemble DESTINATION ${CMAKE_INSTALL_BINDIR})
