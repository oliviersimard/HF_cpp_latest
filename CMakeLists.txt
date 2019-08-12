cmake_minimum_required(VERSION 3.02 FATAL_ERROR)

project(HubbardHF)

# important for the use of tuples
SET(CMAKE_CXX_FLAGS "-std=c++0x")

set(LIBRARIES_CUSTOM fftw3 armadillo json_spirit)

set(SOURCES_COMMON
        src/fft.cpp src/susceptibility_utils.cpp
        src/green_utils.cpp src/json_utils.cpp
        src/integral_utils.cpp)

#set(SOURCES_HF hartree2DSure.cpp ${SOURCES_COMMON})
set(SOURCES_MAIN main.cpp ${SOURCES_COMMON})

#------Executables
#add_executable(HF ${SOURCES_HF})
#target_link_libraries(HF ${LIBRARIES_CUSTOM})
add_executable(MAIN ${SOURCES_MAIN})
target_link_libraries(MAIN ${LIBRARIES_CUSTOM})


# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif()

# Now simply link against gtest or gtest_main as needed. Eg
enable_testing()
add_executable(hubbard_test test/Hubbard_unit_tests.cpp ${SOURCES_COMMON})
target_link_libraries(hubbard_test ${GTEST_LIBRARIES} ${LIBRARIES_CUSTOM} pthread)
add_test(NAME hubbard_test COMMAND ./hubbard_test) ## Allows one to enter "make test" and launch the tests!
