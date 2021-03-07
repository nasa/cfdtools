add_executable(polar_interp
  polar_interp.f90
)
target_link_libraries(polar_interp PRIVATE
  cfdtools_global
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_searchlib
)
add_executable(reorder_spokes reorder_spokes.f90)
target_link_libraries(reorder_spokes PRIVATE
  cfdtools_global
  cfdtools_table_io
)
install(FILES README
        DESTINATION doc
        RENAME polar_interp)
install(TARGETS polar_interp reorder_spokes
        DESTINATION ${CMAKE_INSTALL_BINDIR})
