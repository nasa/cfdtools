add_executable(filter_rows filter_rows.f90)
target_link_libraries(filter_rows
  PRIVATE
    cfdtools_global
    cfdtools_progtools
    cfdtools_table_io
)
install(TARGETS filter_rows DESTINATION ${CMAKE_INSTALL_BINDIR})

# $F90 $FFLAGS -o filter_rows \
#    $ROOT/progtools/string_justify.f90 \
#    $ROOT/table_io/table_io.f90 \
#    filter_rows.f90 \
#    $ROOT/progtools/libprmodules.a

