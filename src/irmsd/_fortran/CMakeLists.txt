cmake_minimum_required(VERSION 3.22)
project(irmsd_fortran LANGUAGES C Fortran)

# Add srcs
set(srcs)
add_subdirectory("fmods")
add_subdirectory("crest_mods")

# Build shared library with Fortran kernels
add_library(irmsd_fortran SHARED
    ${srcs}
)

# OpenMP (parallel iRMSD sorting): optional for source builds, which fall back
# to a serial library. Release wheels set IRMSD_REQUIRE_OPENMP=ON (see
# [tool.cibuildwheel] in pyproject.toml) so they cannot be built serial by accident.
option(IRMSD_USE_OPENMP "Build with OpenMP if the Fortran compiler supports it" ON)
option(IRMSD_REQUIRE_OPENMP "Stop with an error if OpenMP is not found" OFF)
if (IRMSD_REQUIRE_OPENMP AND NOT IRMSD_USE_OPENMP)
    message(FATAL_ERROR "IRMSD_REQUIRE_OPENMP=ON needs IRMSD_USE_OPENMP=ON")
endif()
if (IRMSD_USE_OPENMP)
    if (IRMSD_REQUIRE_OPENMP)
        find_package(OpenMP REQUIRED COMPONENTS Fortran)
    else()
        find_package(OpenMP COMPONENTS Fortran)
    endif()
endif()
if (OpenMP_Fortran_FOUND)
    target_link_libraries(irmsd_fortran PRIVATE OpenMP::OpenMP_Fortran)
    message(STATUS "irmsd: building with OpenMP ${OpenMP_Fortran_VERSION}")
elseif (NOT IRMSD_REQUIRE_OPENMP)
    message(STATUS "irmsd: building without OpenMP (serial sorting)")
endif()

set_target_properties(irmsd_fortran PROPERTIES
    OUTPUT_NAME "irmsd_fortran"
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    POSITION_INDEPENDENT_CODE ON
)

# Install the shared library into the Python package folder inside the wheel
install(TARGETS irmsd_fortran
    LIBRARY DESTINATION irmsd
    RUNTIME DESTINATION irmsd
    ARCHIVE DESTINATION irmsd
)

