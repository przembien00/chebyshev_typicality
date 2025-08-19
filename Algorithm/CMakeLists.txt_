cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 20)


# things to be adapted
set(USE_TYPE DOUBLE)  # specify the type to be used (FLOAT or DOUBLE)
# include_directories(/data/bieniek)  # path to blaze
include_directories(~/Projects/libs/lambda-lanczos/include/) # path to lanczos
#add_definitions(-DEIGEN)  # remove #, if eigen instead of lapack should be used for diagonalization


# Project name
project(chebyshev_typicality)

# Set the default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()

# Find MPI
find_package(MPI REQUIRED)
if(MPI_FOUND)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    include_directories(${MPI_INCLUDE_PATH})
    set(LINK_LIBS ${LINK_LIBS} ${MPI_LIBRARIES})
endif()

# Find HDF5
find_package(HDF5 REQUIRED)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
endif()

# Find Boost
find_package(Boost REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    set(LINK_LIBS ${LINK_LIBS} ${Boost_LIBRARIES})
endif()

# Find LAPACK
# find_package(LAPACK REQUIRED)
# if(LAPACK_FOUND)
#     set(LINK_LIBS ${LINK_LIBS} ${LAPACK_LIBRARIES})
# endif()

# Add definitions
add_definitions(-DUSE_${USE_TYPE})

# Source files
set(SOURCES
    main.cpp
    Functions/Functions.cpp
    Parameter_Space/Parameter_Space.cpp
    Hamiltonians/Hamiltonians.cpp
    Storage_Concept/Storage_Concept.cpp
)


# Header files (for dependencies)
set(GLOB_RECURSE HEADERS
    "main_header.h"
    "Types/*.h"
    "Parameter_Space/*.h"
    "Functions/*.h"
    "Storage_Concept/*.h"
    "cpp_libs/*.h"
)

# Add executable
add_executable(chebyshev_${USE_TYPE} ${SOURCES} ${HEADERS})

# Set the output directory for the executable
set_target_properties(chebyshev_${USE_TYPE} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../
    OUTPUT_NAME ../executable_${USE_TYPE}.out
)

# Link libraries
target_link_libraries(chebyshev_${USE_TYPE} ${LINK_LIBS})