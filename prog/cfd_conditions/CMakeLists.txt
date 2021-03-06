add_executable(cfd_conditions
    earth_atmosphere.f90
    vel_alt_from_Mach_Q.f90
    cfd_conditions.f90
)
target_link_libraries(cfd_conditions PRIVATE
    cfdtools_global
    optlib
    numodules
    progtools
    table_io
)
install(TARGETS cfd_conditions
        DESTINATION ${CMAKE_INSTALL_BINDIR})

