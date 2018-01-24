#-- require at least cmake version 3.0
cmake_minimum_required (VERSION 3.0)

#-- name the project and set the language
project (SuperNu)

#-- set language and find MPI
enable_language (Fortran)
#find_package (MPI REQUIRED) # require MPI for now

#-- set the project version number
set (SuperNu_VERSION_MAJOR 3)
set (SuperNu_VERSION_MINOR 0)

#-- copy dummy commit version file, for now
configure_file(version_dummy.inc version.inc COPYONLY)

#-- add the build directory to the include-file search path
include_directories(${PROJECT_BINARY_DIR})

#-- set the Fortran module directory
set (CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

#-- create lists of top-level source files
file (GLOB_RECURSE SuperNu_FFILES "*.f")
file (GLOB_RECURSE SuperNu_F90FILES "*.f90")

#-- remove versions of mpimod from .f list
file (GLOB MPIMOD "${PROJECT_SOURCE_DIR}/mpimod_mpi.f" "${PROJECT_SOURCE_DIR}/mpimod_ser.f")
configure_file (mpimod_ser.f mpimod.f COPYONLY)
list (REMOVE_ITEM SuperNu_FFILES ${MPIMOD})
list (APPEND SuperNu_FFILES "${PROJECT_BINARY_DIR}/mpimod.f")

#-- add the executable and link MPI libraries
add_executable (supernu ${SuperNu_FFILES} ${SuperNu_F90FILES})

#-- link library and include header for MPI
target_link_libraries (supernu PRIVATE ${MPI_C_LIBRARIES})
target_include_directories (supernu PRIVATE ${MPI_C_INCLUDE_PATH})
