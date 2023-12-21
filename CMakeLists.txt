# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project( LCC_implicit )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Qt5 )

if(CGAL_Qt5_FOUND)
  create_single_source_cgal_program( "test.cpp" )
  target_link_libraries(test PUBLIC CGAL::CGAL_Basic_viewer)
endif()
