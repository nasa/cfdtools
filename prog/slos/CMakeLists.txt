add_executable(slos slos.f90)
target_link_libraries(slos
  PRIVATE
    global
    adt
    eigenlib
    geomlib
    gridlib
    interp3d
    intrinsics
    linsys
    numodules
    progtools
    optlib
    searchlib
    triangulation_io
    ugridlib
    xyzq_io
)
install(TARGETS slos DESTINATION ${CMAKE_INSTALL_BINDIR})

