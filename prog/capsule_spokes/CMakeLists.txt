add_executable(capsule_spokes capsule_spokes.f90)
target_link_libraries(capsule_spokes
  PRIVATE
    cfdtools_global
    cfdtools_geomlib
    cfdtools_numodules
    cfdtools_progtools
    cfdtools_table_io
)
install(TARGETS capsule_spokes DESTINATION ${CMAKE_INSTALL_BINDIR})

