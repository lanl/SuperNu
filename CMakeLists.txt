#-- require at least cmake version 3.8
cmake_minimum_required (VERSION 3.8)

#-- name the project and set the language
project (SuperNu)

#-- set language and find MPI
enable_language (Fortran)
# find_package (MPI REQUIRED)
configure_file (mpimod_ser.f mpimod.f COPYONLY)

#-- set the project version number
set (SuperNu_VERSION_MAJOR 3)
set (SuperNu_VERSION_MINOR 0)

#-- add the build directory to the include-file search path
include_directories(${PROJECT_BINARY_DIR})

#-- copy dummy commit version file, for now
configure_file(version_dummy.inc version.inc COPYONLY)

#-- set the Fortran module directory
set (CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

#-- create lists of top-level source files
file (GLOB_RECURSE SuperNu_FILES "*.f" "*.f90")

#-- remove versions of mpimod from .f list
file (GLOB MPIMOD "mpimod_mpi.f" "mpimod_ser.f")
list (REMOVE_ITEM SuperNu_FILES ${MPIMOD})
list (APPEND SuperNu_FILES "${PROJECT_BINARY_DIR}/mpimod.f")

#-- add the executable and link MPI libraries
add_executable (supernu ${SuperNu_FILES})

#-- link library and include header for MPI
target_link_libraries (supernu PRIVATE ${MPI_C_LIBRARIES})
target_include_directories (supernu PRIVATE ${MPI_C_INCLUDE_PATH})
