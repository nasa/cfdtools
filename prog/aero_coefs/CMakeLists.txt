add_executable(aero_coefs
    aero_coefs.f90
    patch_coefs.f90
)

target_link_libraries(aero_coefs
  PRIVATE
    global
    gridlib
    geomlib
    progtools
    xyzq_io
)

install(TARGETS aero_coefs DESTINATION ${CMAKE_INSTALL_BINDIR})

