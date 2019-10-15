add_executable(free_stream_conditions free_stream_conditions.f90)
target_link_libraries(free_stream_conditions
  PRIVATE
    global
    progtools
    table_io
)
install(TARGETS free_stream_conditions DESTINATION ${CMAKE_INSTALL_BINDIR})

# $F90 $FFLAGS -o free_stream_conditions \
#    $ROOT/progtools/string_justify.f90 \
#    $ROOT/table_io/table_io.f90 \
#    free_stream_conditions \
#    $ROOT/progtools/libprmodules.a

