#-- require at least cmake version 3.8
cmake_minimum_required (VERSION 3.8)

#-- specify configuration options
option (USE_MPI "USE_MPI" ON)
option (USE_OPENMP "USE_OPENMP" ON)
option (USE_HDF5 "USE_HDF5" OFF)

#-- name the project and set the language
project (SuperNu)
set (SuperNu_VERSION_MAJOR 4)
set (SuperNu_VERSION_MINOR "x")
enable_language (Fortran)

#-- add the build directory to the include-file search path
include_directories(${PROJECT_BINARY_DIR})

#-- copy dummy commit version file, for now
configure_file(version_dummy.inc version.inc COPYONLY)

#-- find MPI (cmake -DUSE_MPI=ON <src dir>)
find_package (MPI)
if (MPI_Fortran_FOUND AND USE_MPI)
  configure_file (mpimod_mpi.f mpimod.f COPYONLY)
else ()
  configure_file (mpimod_ser.f mpimod.f COPYONLY)
endif ()

#-- find OpenMP (cmake -DUSE_OPENMP=ON <src dir>)
find_package (OpenMP)
if (OPENMP_FOUND AND USE_OPENMP)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

#-- find HDF5 with Fortran bindings
if (USE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS Fortran)
endif ()

#-- set the Fortran module directory
set (CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

#-- create lists of top-level source files
file (GLOB_RECURSE SuperNu_FILES "*.f" "*.f90")
if (NOT USE_HDF5)
  list (REMOVE_ITEM SuperNu_FILES "${PROJECT_SOURCE_DIR}/h5aux.f90")
endif ()

#-- remove versions of mpimod from .f list
file (GLOB MPIMOD "mpimod_mpi.f" "mpimod_ser.f")
list (REMOVE_ITEM SuperNu_FILES ${MPIMOD})
list (APPEND SuperNu_FILES "${PROJECT_BINARY_DIR}/mpimod.f")

#-- add the executable and link MPI libraries
add_executable (supernu ${SuperNu_FILES})

#-- link library and include header for MPI
if (MPI_Fortran_FOUND AND USE_MPI)
  target_link_libraries (supernu PRIVATE ${MPI_C_LIBRARIES})
  target_include_directories (supernu PRIVATE ${MPI_C_INCLUDE_PATH})
  target_link_libraries (supernu PRIVATE ${MPI_Fortran_LIBRARIES})
  target_include_directories (supernu PRIVATE ${MPI_Fortran_INCLUDE_PATH})
endif ()

if (HDF5_FOUND AND USE_HDF5)
  target_include_directories(supernu PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(supernu PRIVATE ${HDF5_Fortran_LIBRARIES})
endif ()

#-- indicate successful configuration
message ("")
message ("SuperNu Version " ${SuperNu_VERSION_MAJOR} "."
	${SuperNu_VERSION_MINOR} " configured.")
message ("")
