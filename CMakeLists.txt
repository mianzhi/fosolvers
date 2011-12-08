CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(FOSolverS)
ENABLE_LANGUAGE(Fortran)

# compiler
SET(CMAKE_Fortran_COMPILER mpif90)
ADD_DEFINITIONS(-Wall -O3)
SET(CMAKE_EXE_LINKER_FLAGS -Wall -fopenmp)
SET(CMAKE_SHARED_LINKER_FLAGS -Wall -fopenmp)

# additional include/library directories
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/libfosolvers)
LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/libfosolvers)

# subdirectories
ADD_SUBDIRECTORY(libfosolvers)
