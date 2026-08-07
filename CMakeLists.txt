# * *************************************************************
# *  
# *   Soft Active Mater on Surfaces (SAMoS)
# *   
# *   Author: Rastko Sknepnek
# *  
# *   Division of Physics
# *   School of Engineering, Physics and Mathematics
# *   University of Dundee
# *   
# *   (c) 2013, 2014
# * 
# *   School of Science and Engineering
# *   School of Life Sciences 
# *   University of Dundee
# * 
# *   (c) 2015
# * 
# *   Author: Silke Henkes
# * 
# *   Department of Physics 
# *   Institute for Complex Systems and Mathematical Biology
# *   University of Aberdeen  
# * 
# *   (c) 2014, 2015
# *  
# *   This program cannot be used, copied, or modified without
# *   explicit written permission of the authors.
# * 
# * ************************************************************** 

cmake_minimum_required(VERSION 3.16...3.30)

# Policies needed by older dependency configuration files
if(POLICY CMP0003)
    cmake_policy(SET CMP0003 NEW)
endif()

# Temporary compatibility for dependency packages that still use FindBoost.
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

project(SAMoS LANGUAGES C CXX)

# CGAL requires C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ENABLE_STATIC)
  SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  SET (THREAD_LIB -lpthread)
endif(ENABLE_STATIC)

# Setup a number of misc options and libraries
include (CMakeMiscSetup.txt)
# Find the boost libraries and set them up
include (CMakeBoostSetup.txt)
# Find GSL libraries
include (CMakeGSLSetup.txt)
# Set default CFlags
include (CMakeCFlagsSetup.txt)
# Configure some source files, include directories, and create variables listing all source files
include (CMakeSRCSetup.txt)
# Configure VTK libraries
include (CMakeVTKSetup.txt)
# Configure CGAL libraries
include (CMakeCGALSetup.txt)

################################
## Define common libraries used by every target in MEMBRANE
set(BOOST_LIBS 	${Boost_REGEX_LIBRARY} ${Boost_IOSTREAMS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} )
set(GSL_LIBS ${GSL_LIBRARIES})
set(MATH_LIB -lm)
set(SAMoS_LIBS ${GSL_LIBS} ${MATH_LIB} ${BOOST_LIBS} ${VTK_LIBS})

# ##############################################
# place all executables in the build directory 
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

###############################
## include documentation directories
if (ENABLE_DOXYGEN)
	add_subdirectory (doc)
endif (ENABLE_DOXYGEN)

add_subdirectory(src)

