CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

PROJECT(Clotting LANGUAGES C CXX)
FIND_PACKAGE(IBAMR REQUIRED)
FIND_PACKAGE(ADS REQUIRED)

SET(UNBOUND_CXX_SRC
    src/BondUnactivatedSource.cpp
    src/BoundaryMeshMapping.cpp
    src/CohesionStressUnactivatedRHS.cpp
    src/UnactivatedPlateletSource.cpp
    src/utility_functions.cpp
    )

SET(BOUND_CXX_SRC
    src/BondBoundSource.cpp
    src/BoundPlateletSource.cpp
    src/BoundaryMeshMapping.cpp
    src/BoundVelocitySource.cpp
    src/CellToFaceFcn.cpp
    src/ClotParameters.cpp
    src/CohesionStressBoundRHS.cpp
    src/DragForce.cpp
    src/utility_functions.cpp
    )

ADD_LIBRARY(UnboundClot)
TARGET_SOURCES(UnboundClot PRIVATE ${UNBOUND_CXX_SRC})
TARGET_INCLUDE_DIRECTORIES(UnboundClot PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:${CMAKE_INSTALL_DIR}/include/>)
TARGET_LINK_LIBRARIES(UnboundClot PUBLIC IBAMR::IBAMR2d)
TARGET_LINK_LIBRARIES(UnboundClot PUBLIC ADS::ADS2d)

ADD_LIBRARY(BoundClot)
TARGET_SOURCES(BoundClot PRIVATE ${BOUND_CXX_SRC})
TARGET_INCLUDE_DIRECTORIES(BoundClot PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:${CMAKE_INSTALL_DIR}/include/>)
TARGET_LINK_LIBRARIES(BoundClot PUBLIC IBAMR::IBAMR2d)
TARGET_LINK_LIBRARIES(BoundClot PUBLIC ADS::ADS2d)

ADD_SUBDIRECTORY(fpc)
ADD_SUBDIRECTORY(aneurysm)
ADD_SUBDIRECTORY(tests)
ADD_SUBDIRECTORY(channel)
