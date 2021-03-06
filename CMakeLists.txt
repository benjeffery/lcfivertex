########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, DESY
########################################################


########################################################
# CMake compatibility issues: don't modify this, please!
CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
########################################################



####################################################################################################
# User section: modify as needed
####################################################################################################


# project name
PROJECT( LCFIVertex )

# project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 2 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 2 )


### SETTINGS #################################################################

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )

# project dependencies
# e.g. SET( ${PROJECT_NAME}_DEPENDS "Marlin MarlinUtil LCIO GEAR CLHEP GSL RAIDA" )
SET( ${PROJECT_NAME}_DEPENDS "Marlin MarlinUtil LCIO GEAR" )

# set default cmake build type to RelWithDebInfo
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

# set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()


### INCLUDES #################################################################

# include directories
INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/include" )
INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/vertex_lcfi" )
INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/vertex_lcfi/nnet/inc" )
INCLUDE_DIRECTORIES( "diagnostics/include" )
# boost
INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/boost" )

# install include files
INSTALL( DIRECTORY
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/boost"
        DESTINATION .
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )
INSTALL( DIRECTORY
        "${PROJECT_SOURCE_DIR}/diagnostics/include"
        DESTINATION "diagnostics"
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )
INSTALL( DIRECTORY
        "${PROJECT_SOURCE_DIR}/vertex_lcfi/inc"
        DESTINATION "vertex_lcfi"
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )
INSTALL( DIRECTORY
        "${PROJECT_SOURCE_DIR}/vertex_lcfi/algo/inc"
        DESTINATION "vertex_lcfi/algo"
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )
INSTALL( DIRECTORY
        "${PROJECT_SOURCE_DIR}/vertex_lcfi/nnet/inc"
        DESTINATION "vertex_lcfi/nnet"
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )
INSTALL( DIRECTORY
        "${PROJECT_SOURCE_DIR}/vertex_lcfi/util/inc"
        DESTINATION "vertex_lcfi/util"
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )
INSTALL( DIRECTORY
        "${PROJECT_SOURCE_DIR}/vertex_lcfi/zvtop/include"
        DESTINATION "vertex_lcfi/zvtop"
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE )


### SOURCES ##################################################################

# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )

# add debug definitions
#IF( CMAKE_BUILD_TYPE STREQUAL "Debug" OR
#    CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
#    ADD_DEFINITIONS( "-DDEBUG" )
#ENDIF()

# get list of all source files
AUX_SOURCE_DIRECTORY( src processor_srcs )
AUX_SOURCE_DIRECTORY( diagnostics/src diagnostics_srcs )
FILE( GLOB_RECURSE vertex_lcfi_srcs "vertex_lcfi/*.cpp" )

# library sources
SET( library_sources ${processor_srcs} ${diagnostics_srcs} ${vertex_lcfi_srcs} )


### DOCUMENTATION ############################################################

FIND_PACKAGE( Doxygen )
FIND_PACKAGE( LATEX )

# global doc target for calling doc_html and doc_latex
ADD_CUSTOM_TARGET( doc )

#doxygen
IF( DOXYGEN_FOUND )
    
    ADD_CUSTOM_COMMAND(
        OUTPUT  "${PROJECT_SOURCE_DIR}/doc/html"
        COMMAND "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
        COMMENT "Building API Documentation..."
        VERBATIM )
    
    # add doc_html target
    ADD_CUSTOM_TARGET( doc_html DEPENDS
        "${PROJECT_SOURCE_DIR}/doc/html" )
    # tell doc target to call doc_html
    ADD_DEPENDENCIES( doc doc_html )
    
ELSE()
    MESSAGE( STATUS "Doxygen not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()

#latex
IF( DOXYGEN_FOUND AND LATEX_COMPILER AND
    MAKEINDEX_COMPILER AND PDFLATEX_COMPILER )
     
     # generate latex with doxygen
     ADD_CUSTOM_COMMAND(
        OUTPUT  "${PROJECT_SOURCE_DIR}/doc/latex"
        COMMAND "${DOXYGEN_EXECUTABLE}" ARGS Doxyfile.latex
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
        COMMENT "Generating Latex Documentation with Doxygen..."
        VERBATIM )
    
    # run latex + pdflatex
    ADD_CUSTOM_COMMAND(
        OUTPUT  "${PROJECT_SOURCE_DIR}/doc/latex/refman.pdf"
        COMMAND "${CMAKE_BUILD_TOOL}"
        # copy pdf to doc directory
        COMMAND "${CMAKE_COMMAND}" -E copy "refman.pdf" "${PROJECT_SOURCE_DIR}/doc"
        DEPENDS "${PROJECT_SOURCE_DIR}/doc/latex"
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc/latex"
        COMMENT "Building Latex Documentation..."
        VERBATIM )
    
    # add doc_latex target
    ADD_CUSTOM_TARGET( doc_latex DEPENDS
        "${PROJECT_SOURCE_DIR}/doc/latex/refman.pdf" )
    # tell doc target to call doc_latex
    ADD_DEPENDENCIES( doc doc_latex )
    
ELSE()
    MESSAGE( STATUS "Latex not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()

# install documentation
IF( INSTALL_DOC )
    # make sure doxygen is executed (make doc) before make install
    INSTALL( CODE "EXEC_PROGRAM(${CMAKE_BUILD_TOOL} ${PROJECT_BINARY_DIR} ARGS doc)" )

    # install documentation
    INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
            DESTINATION .
            PATTERN "*CVS*" EXCLUDE )
ENDIF()


##########################################################################################
# End of User section: please try not to modify below this line
##########################################################################################


# library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION
    "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )
SET( ${PROJECT_NAME}_VERSION
    "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )

# add install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# add install path to the rpath list (apple)
IF( APPLE )
    SET( CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
    MARK_AS_ADVANCED( CMAKE_INSTALL_NAME_DIR )
ENDIF()

# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

# output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# DEPENDENCIES: this code has to be placed before adding any library or
# executable so that these are linked properly against the dependencies
IF( DEFINED ${PROJECT_NAME}_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
    IF( NOT DEFINED CMAKE_MODULE_PATH )
        MESSAGE( FATAL_ERROR "CMAKE_MODULE_PATH not set! Set it with: "
            "-DCMAKE_MODULE_PATH=\"/path_to_ilcsoft/CMakeModules\"" )
    ENDIF()
    # load macro
    INCLUDE( "MacroCheckDeps" )
    CHECK_DEPS()
ENDIF()

# LIBRARY
ADD_LIBRARY( lib_${PROJECT_NAME} ${library_sources} )
# create symbolic lib target for calling target lib_XXX
ADD_CUSTOM_TARGET( lib DEPENDS lib_${PROJECT_NAME} )
# change lib_target properties
SET_TARGET_PROPERTIES( lib_${PROJECT_NAME} PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${${PROJECT_NAME}_VERSION}
    SOVERSION ${${PROJECT_NAME}_SOVERSION}
    # allow creating static and shared libs without conflicts
    CLEAN_DIRECT_OUTPUT 1
    # avoid conflicts between library and binary target names
    OUTPUT_NAME ${PROJECT_NAME} )

# install library
INSTALL( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

# create uninstall configuration file 
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" )

# create configuration file from .in file
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in"
                "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY )

# install configuration file
INSTALL( FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" DESTINATION . )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}" )
MESSAGE( STATUS "${PROJECT_NAME}_DEPENDS = \"${${PROJECT_NAME}_DEPENDS}\"" )
MESSAGE( STATUS "BUILD_WITH = \"${BUILD_WITH}\"" )
MESSAGE( STATUS "INSTALL_DOC = ${INSTALL_DOC}" )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

# force some variables that could be defined in the command line
# to be written to cache
SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
    "Set to OFF to build static libraries" FORCE )
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
    "Where to install ${PROJECT_NAME}" FORCE )
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
SET( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" CACHE PATH
    "Path to custom CMake Modules" FORCE )
SET( INSTALL_DOC "${INSTALL_DOC}" CACHE BOOL
    "Set to OFF to skip build/install Documentation" FORCE )

# export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "${PROJECT_NAME}BuildSettings.cmake" )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "${PROJECT_NAME}LibDeps.cmake" )

