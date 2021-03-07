add_executable(ensemble
  ensemble.f90
)
target_link_libraries(ensemble cfdtools_global)
target_link_libraries(ensemble
    cfdtools_intrinsics
    cfdtools_tecplot_adt
    cfdtools_tecplot_gridlib
    cfdtools_tecplot_io
    cfdtools_tecplot_surface_patch_utilities
    cfdtools_progtools
)
install(TARGETS ensemble DESTINATION ${CMAKE_INSTALL_BINDIR})
