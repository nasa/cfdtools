add_executable(capsule_spokes capsule_spokes.f90)
target_link_libraries(capsule_spokes
  PRIVATE
    global
    geomlib
    numodules
    progtools
    table_io
)
install(TARGETS capsule_spokes DESTINATION ${CMAKE_INSTALL_BINDIR})

