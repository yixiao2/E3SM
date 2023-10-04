# - Try to find PIO
#
# Placeholder code until spio gets a proper config.cmake
#
# Once done, this will define:
#
#   The "spio" target
#

if (TARGET spio)
  return()
endif()

if (NOT PIO_LIBDIR)
  set(PIO_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

# All pio libs will link against netcdf libraries
find_package(NETCDF REQUIRED)

if (PIO_VERSION STREQUAL 2)
  # This is a pio2 library
  set (PIOLIBS pioc piof)

  # Wrap pioc lib with its deps in a target
  add_library (pioc STATIC IMPORTED GLOBAL)
  set_property(TARGET pioc PROPERTY
    IMPORTED_LOCATION "${PIO_LIBDIR}/libpioc.a")
  target_link_libraries(pioc INTERFACE netcdf)

  # Not all machines/PIO installations use ADIOS but, for now,
  # we can assume that an MPI case with ADIOS2_ROOT set is probably
  # using adios.
  if (NOT MPILIB STREQUAL mpi-serial AND DEFINED ENV{ADIOS2_ROOT})
    find_package(MPI REQUIRED COMPONENTS C)
    find_package(ADIOS2 REQUIRED COMPONENTS C)
    target_link_libraries(pioc INTERFACE adios2::adios2)
  endif()

  # Wrap piof lib with its deps in a target
  add_library (piof STATIC IMPORTED GLOBAL)
  set_property(TARGET piof PROPERTY
    IMPORTED_LOCATION "${PIO_LIBDIR}/libpioc.a")
  target_link_libraries(piof INTERFACE pioc netcdf)
else()
  # This is a pio1 library
  set (PIOLIBS pio)

  add_library (pio1 STATIC IMPORTED GLOBAL)
  set_property(TARGET pio1 PROPERTY
    IMPORTED_LOCATION "${PIO_LIBDIR}/libpio.a")
  target_link_libraries(pio INTERFACE netcdf)
endif()

# Declare a single spio target, since that's expected downstream
# Ideally, each downstream app should just use pioc/piof, based
# on their needs, but we may as well provide a "generic" spio target.
add_library(spio INTERFACE)
target_link_libraries(spio INTERFACE ${PIOLIBS})
