cmake_minimum_required(VERSION 3.0)

PROJECT(mic_array)

# headfiles 
INCLUDE_DIRECTORIES(include)


# src
 AUX_SOURCE_DIRECTORY(src DIR_SRCS)

SET(TEST_MATH ${DIR_SRCS})

# LINK_LIBRARIES(m) 

ADD_EXECUTABLE(${PROJECT_NAME} ${TEST_MATH})