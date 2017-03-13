CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(fabili)

SET(CMAKE_CXX_COMPILER "/opt/rh/devtoolset-3/root/usr/bin/g++")
SET(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR} )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -Wextra -Werror -O3")


SET( INCLUDE_DIRS 
    ${INCLUDE_DIRS}
    ${PROJECT_BINARY_DIR}  
    ${PROJECT_SOURCE_DIR}  
    $ENV{EIGEN_INCLUDE_DIR}
)

INCLUDE_DIRECTORIES( ${INCLUDE_DIRS} )

ADD_EXECUTABLE(run.exe run.cxx)

ADD_LIBRARY(fabili      SHARED fabili.cc         )
ADD_LIBRARY(testclass   SHARED testclass.cc      )

TARGET_LINK_LIBRARIES( testclass fabili )

TARGET_LINK_LIBRARIES( run.exe fabili )
TARGET_LINK_LIBRARIES( run.exe testclass )
TARGET_LINK_LIBRARIES( run.exe testclass )
