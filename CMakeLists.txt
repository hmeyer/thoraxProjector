cmake_minimum_required(VERSION 2.6)
PROJECT(thoraxProjector)

FIND_PATH(MYLIBS itkbasics.h ../myLibs ../../myLibs)

SET(SRC ${MYLIBS}/itkbasics.cpp)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)



find_package(Boost REQUIRED filesystem)

SET(Libraries
 ITKCommon
 ITKAlgorithms
 ITKIO
)

SET(CurrentExe "thoraxProjector")
ADD_EXECUTABLE(${CurrentExe} thoraxProj.cpp ${SRC})
INCLUDE_DIRECTORIES(${MYLIBS})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries} ${Boost_LIBRARIES})


