#  CMAKE CONFIG
CMAKE_MINIMUM_REQUIRED( VERSION 3.5 FATAL_ERROR )
project(RapidSim)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(RAPIDSIM_ROOT "${CMAKE_INSTALL_PREFIX}" CACHE STRING
    "Directory to install RapidSim's data to")

#This module respects HFS, e.g. defines lib or lib64 when it is needed.
include("GNUInstallDirs")


SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
SET( ${PROJECT_NAME}_MINOR_VERSION 4 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )

SET( CMAKE_BUILD_TYPE Release )

EXECUTE_PROCESS(
            COMMAND git rev-parse HEAD
        OUTPUT_VARIABLE GIT_REV )
MESSAGE( "GIT VERSION:  ${GIT_REV}" )

EXECUTE_PROCESS(
        COMMAND date +%Y_%m_%d
        OUTPUT_VARIABLE BUILD_DATE )

SET ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} )

if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "i686")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -m3dnow")
endif()

# Add the other flags regardless of architecture
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fmerge-all-constants -D__ROOFIT_NOBANNER -Wall -Wextra -Werror")

set(CMAKE_CXX_FLAGS_DEBUG          "-g3 -Wall")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g3 -O2")


# setting linker flags
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(CMAKE_EXE_LINKER_FLAGS "-Wl,-dead_strip_dylibs")
  SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,-dead_strip_dylibs")
  SET(STATIC_LIBRARY_FLAGS "-Wl,-dead_strip_dylibs")
ELSEIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed")
  SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--as-needed")
  SET(STATIC_LIBRARY_FLAGS "-Wl,--as-needed")
ENDIF()

find_package(ROOT CONFIG REQUIRED COMPONENTS Core RIO RooFit RooFitCore RooStats Hist Tree Matrix Physics MathCore)
include(${ROOT_USE_FILE})
message(STATUS "ROOT includes: ${ROOT_INCLUDE_DIRS}")
message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")
message(STATUS "ROOT library directory: ${ROOT_LIBRARY_DIR}")
include_directories(${ROOT_INCLUDE_DIRS})

find_package(EvtGen CONFIG)
if(EvtGen_FOUND)
  message(STATUS "Found EvtGen: ${EvtGen_DIR}")
  message(STATUS "Will attempt to link against EvtGen")
  set(EVTGEN    EvtGen::EvtGen)
  set(EVTGENEXT EvtGen::EvtGenExternal)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRAPID_EVTGEN")
  if (CMAKE_CXX_COMPILER_ID MATCHES Clang AND CMAKE_SYSTEM_NAME MATCHES Darwin)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-undefined,dynamic_lookup")
  else()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--unresolved-symbols=ignore-in-shared-libs")
  endif()
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
elseif(DEFINED ENV{EVTGEN_ROOT})
  message(STATUS "EVTGEN_ROOT is set: ${EVTGEN_ROOT}")
  message(STATUS "Will attempt to link against EvtGen")
  find_library(EVTGEN EvtGen PATHS "$ENV{EVTGEN_ROOT}/lib")
  find_library(EVTGENEXT EvtGenExternal PATHS "$ENV{EVTGEN_ROOT}/lib")
  if(NOT EVTGEN)
    message(INFO " Will not link against EvtGen : libEvtGen not found")
  else()
    if(NOT EVTGENEXT)
      message(INFO " Will not link against EvtGen : libEvtGenExternal not found")
    else()
      include_directories(${EVTGEN_ROOT})
      set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} -DRAPID_EVTGEN -I$ENV{EVTGEN_ROOT}")
      if (CMAKE_CXX_COMPILER_ID MATCHES Clang AND CMAKE_SYSTEM_NAME MATCHES Darwin)
        # Apple's build of lld doesn't support --unresolved-symbols
        set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
      else()
        set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} --unresolved-symbols=ignore-in-shared-libs")
      endif()
      set(EvtGen_FOUND TRUE)
    endif()
  endif()
else()
  message(INFO " Will not link against EvtGen : EVTGEN_ROOT not defined")
endif()

add_subdirectory(src)

set(bindir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

if (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  configure_file(${CMAKE_SOURCE_DIR}/validation/runValidation.sh.in ${RAPIDSIM_ROOT}/${CMAKE_INSTALL_BINDIR}/runValidation.sh   @ONLY)
endif()

install(DIRECTORY validation DESTINATION ${RAPIDSIM_ROOT} )
install(DIRECTORY rootfiles DESTINATION ${RAPIDSIM_ROOT} )
install(DIRECTORY utils DESTINATION ${RAPIDSIM_ROOT} )
install(DIRECTORY config DESTINATION ${RAPIDSIM_ROOT} )
