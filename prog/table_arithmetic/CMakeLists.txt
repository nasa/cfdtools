add_executable(table_arithmetic table_arithmetic.f90)
target_link_libraries(table_arithmetic
  PRIVATE
    cfdtools_global
    cfdtools_numodules
    cfdtools_progtools
    cfdtools_searchlib
    cfdtools_table_io
)
install(TARGETS table_arithmetic DESTINATION ${CMAKE_INSTALL_BINDIR})


# $F90 $FFLAGS -o table_arithmetic \
#    $ROOT/progtools/string_justify.f90 \
#    $ROOT/table_io/table_io.f90 \
#    table_arithmetic.f90 \
#    $ROOT/numodules/f90numodules.a \
#    $ROOT/searchlib/f90searchlib.a \
#    $ROOT/progtools/libprmodules.a

