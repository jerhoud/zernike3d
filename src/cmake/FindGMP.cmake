include(FindPackageHandleStandardArgs)

find_library(GMP_C_LIBRARIES
  NAMES libgmp.a gmp
  DOC "GMP C libraries"
  PATHS $ENV{GMPDIR}
)
find_library(GMP_CXX_LIBRARIES
  NAMES libgmpxx.a gmpxx
  DOC "GMP C++ libraries"
  PATHS $ENV{GMPDIR}
)

find_path(GMP_C_INCLUDES
  NAMES gmp.h
  DOC "GMP C header"
  PATHS $ENV{GMPDIR}
)

find_path(GMP_CXX_INCLUDES
  NAMES gmpxx.h
  DOC "GMP C++ header"
  PATHS $ENV{GMPDIR}
)

find_package_handle_standard_args(GMP
	REQUIRED_VARS GMP_C_LIBRARIES GMP_C_INCLUDES GMP_CXX_LIBRARIES GMP_CXX_INCLUDES)

if (GMP_FOUND)
  if (NOT TARGET GMP::GMP)
    add_library(GMP::GMP UNKNOWN IMPORTED)
    set_target_properties(GMP::GMP PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${GMP_C_INCLUDES}"
      IMPORTED_LOCATION "${GMP_C_LIBRARIES}")
  endif()

  if (NOT TARGET GMP::GMPXX)
    add_library(GMP::GMPXX UNKNOWN IMPORTED)
    set_target_properties(GMP::GMPXX PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${GMP_CXX_INCLUDES}"
      IMPORTED_LOCATION "${GMP_CXX_LIBRARIES}")
    target_link_libraries(GMP::GMPXX INTERFACE GMP::GMP)
  endif()
endif()
