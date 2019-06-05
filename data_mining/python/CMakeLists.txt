#python modules
find_package(Boost 1.67.0 COMPONENTS python36 timer REQUIRED)
find_package(PythonLibs 3.6 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
message(STATUS "Boost_VERSION: ${Boost_VERSION}")
message(STATUS "Python_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Python_LIBRARIES: ${PYTHON__LIBRARY}")

include_directories(${CMAKE_SOURCE_DIR})

# Build our library
add_library(graphgenerator SHARED ${CMAKE_SOURCE_DIR}/GraphGenerator.cxx)

# Define the wrapper library that wraps our library
add_library(tacos SHARED
        GraphGenerator_Python.cxx
        ${CMAKE_SOURCE_DIR}/MainTable.cxx
        ${CMAKE_SOURCE_DIR}/SamplingModel.cxx
        ${CMAKE_SOURCE_DIR}/RandomGen.cxx
        ${CMAKE_SOURCE_DIR}/LabelsReshuffler.cxx)
# don't prepend wrapper library name with lib
set_target_properties(tacos PROPERTIES PREFIX "" )
# must be .so (not .dylib)
set_target_properties(tacos PROPERTIES SUFFIX ".so")

target_link_libraries(tacos ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} graphgenerator curl)
