CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

PROJECT(Clotting LANGUAGES C CXX)
FIND_PACKAGE(IBAMR REQUIRED)
FIND_PACKAGE(ADS REQUIRED)

SET(CXX_SRC
    src/BondSource.cpp
    src/BoundPlateletSource.cpp
    src/BoundaryMeshMapping.cpp
    src/BoundVelocityFunction.cpp
    src/CohesionStressBoundRHS.cpp
    src/CohesionStressUnactivatedRHS.cpp
    src/UnactivatedPlateletSource.cpp
    src/utility_functions.cpp
    )

ADD_LIBRARY(Clot)
TARGET_SOURCES(Clot PRIVATE ${CXX_SRC})
TARGET_INCLUDE_DIRECTORIES(Clot PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:${CMAKE_INSTALL_DIR}/include/>)
TARGET_LINK_LIBRARIES(Clot PUBLIC IBAMR::IBAMR2d)
TARGET_LINK_LIBRARIES(Clot PUBLIC ADS::ADS2d)


ADD_SUBDIRECTORY(fpc)
ADD_SUBDIRECTORY(aneurysm)
ADD_SUBDIRECTORY(tests)
