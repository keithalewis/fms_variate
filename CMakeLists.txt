cmake_minimum_required(VERSION 3.15)

project(fms_variate
	VERSION 0.1
	DESCRIPTION "Random variable distributions"
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
#find_package(GSL REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)

include(CTest)
enable_testing()

add_executable(fms_variate.t)
target_sources(fms_variate.t
	PRIVATE
	fms_variate.t.cpp
	fms_variate_constant.t.cpp
	fms_variate_normal.t.cpp)

add_test(NAME fms_variate.t COMMAND fms_variate.t)
