CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(FOSolverS)
ENABLE_LANGUAGE(Fortran)

# compiler & linker arguments
SET(CMAKE_Fortran_COMPILER mpif90)
SET(CMAKE_Fortran_FLAGS "-Wall -Ofast -fopenmp")
SET(CMAKE_EXE_LINKER_FLAGS "-Wall -Ofast -fopenmp")
SET(CMAKE_SHARED_LINKER_FLAGS "-Wall -Ofast -fopenmp")

# module, executable & lib output directories
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# include/library directories
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/include)
LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/lib)

# subdirectories
ADD_SUBDIRECTORY(src)
