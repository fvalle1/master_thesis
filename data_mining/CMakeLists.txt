cmake_minimum_required(VERSION 3.12)
project(thesis_data_mining)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(BOOST_ROOT ~/boost)

find_package(Boost 1.67.0 COMPONENTS system timer random REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(
            ~/anaconda3/include
            ~/Developer/pcg
    )

    add_executable(
            data_mining
            TelegramWatch.cxx
            main.cpp
            MainTable.cxx
            SamplingModel.cxx
            RandomGen.cxx
            GraphGenerator.cxx
            BioParameters.cxx
    )

    target_link_libraries(data_mining ${Boost_LIBRARIES})
    target_link_libraries(data_mining curl)
endif()

#python modules
find_package(Boost 1.67.0 COMPONENTS python36 timer REQUIRED)
find_package(PythonLibs 3.6 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
message(STATUS "Boost_VERSION: ${Boost_VERSION}")
message(STATUS "Python_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Python_LIBRARIES: ${PYTHON__LIBRARY}")

# Build our library
add_library(graphgenerator SHARED GraphGenerator.cxx)

# Define the wrapper library that wraps our library
add_library(tacos SHARED
        GraphGenerator_Python.cxx
        MainTable.cxx
        SamplingModel.cxx
        RandomGen.cxx)
# don't prepend wrapper library name with lib
set_target_properties(tacos PROPERTIES PREFIX "" )
# must be .so (not .dylib)
set_target_properties(tacos PROPERTIES SUFFIX ".so")

target_link_libraries(tacos ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} graphgenerator curl)
