add_library(adt
  adt_utilities.f90
  adt_utilities.f90.before-structured-curve-case
  new_adt_utilities.f90.unnecessary-rotation
)
target_link_libraries(adt global)
install(FILES README DESTINATION doc)
