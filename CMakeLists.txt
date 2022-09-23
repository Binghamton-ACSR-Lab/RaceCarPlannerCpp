cmake_minimum_required(VERSION 3.16)
project(RaceCarPlanner)

set(CMAKE_CXX_STANDARD 17)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
#add_subdirectory(test)
#enable_testing()

find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)
find_package(SQLiteCpp REQUIRED)
find_package(yaml-cpp)
find_package(Threads)
find_package(TBB REQUIRED)

include_directories(${PYTHON3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})

#set(UNIT_TESTS test/csv_reader_test.cpp)

include_directories(include
        include/wheel
        include/utility
        include/chassis
        include/track
        include/global_planner
        include/vehicle
        include/map)

add_executable(RaceCarPlanner main.cpp)
target_include_directories(RaceCarPlanner PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS} )
target_link_libraries(RaceCarPlanner SQLiteCpp sqlite3 casadi gtest gtest_main ${PYTHON_LIBRARIES} Python3::NumPy yaml-cpp pthread dl tbb tinyxml2)
#add_test( RaceCarPlanner RaceCarPlanner )



