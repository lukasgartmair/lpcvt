cmake_minimum_required (VERSION 2.4) 
project (LpCVT)


IF(WIN32)
	find_package(CGAL QUIET REQUIRED COMPONENTS Core )
ELSE(WIN32)
	find_package(CGAL QUIET COMPONENTS Core )
ENDIF(WIN32)

#
# Find the CppUnit includes and library
#
# This module defines
# CPPUNIT_INCLUDE_DIR, where to find tiff.h, etc.
# CPPUNIT_LIBRARIES, the libraries to link against to use CppUnit.
# CPPUNIT_FOUND, If false, do not try to use CppUnit.

# also defined, but not for general use are
# CPPUNIT_LIBRARY, where to find the CppUnit library.
# CPPUNIT_DEBUG_LIBRARY, where to find the CppUnit library in debug

FIND_PATH(CPPUNIT_INCLUDE_DIR cppunit/TestCase.h
  /usr/local/include
  /usr/include
  /usr/include/cppunit
)

# With Win32, important to have both
IF(WIN32)
  FIND_LIBRARY(CPPUNIT_LIBRARY cppunit
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
  FIND_LIBRARY(CPPUNIT_DEBUG_LIBRARY cppunitd
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
ELSE(WIN32)
  # On unix system, debug and release have the same name
  FIND_LIBRARY(CPPUNIT_LIBRARY cppunit
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
  FIND_LIBRARY(CPPUNIT_DEBUG_LIBRARY cppunit
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
ENDIF(WIN32)

IF(CPPUNIT_INCLUDE_DIR)
  IF(CPPUNIT_LIBRARY)
    SET(CPPUNIT_FOUND "YES")
    SET(CPPUNIT_LIBRARIES ${CPPUNIT_LIBRARY} ${CMAKE_DL_LIBS})
    SET(CPPUNIT_DEBUG_LIBRARIES ${CPPUNIT_DEBUG_LIBRARY}
${CMAKE_DL_LIBS})
    SET(CPPUNIT_DIR /usr/include/cppunit)
    SET(cppunit_DIR /usr/include/cppunit)
  ENDIF(CPPUNIT_LIBRARY)
ENDIF(CPPUNIT_INCLUDE_DIR)

	
IF(CGAL_FOUND)
	MESSAGE("CGAL found   ${CGAL_USE_FILE}")
	INCLUDE( ${CGAL_USE_FILE} )
	SET( CMAKE_VERBOSE_MAKEFILE OFF )
ELSE(CGAL_FOUND)
	IF(WIN32)
		MESSAGE("Error : CGAL not found")
	ELSE(WIN32)
		SET(CGAL_DIR /opt/CGAL-3.4 )
		SET(CGAL_LIBS_DEBUG   CGAL_Core CGAL gmp mpfr boost_thread-mt )	
		SET(CGAL_LIBS_RELEASE CGAL_Core CGAL gmp mpfr boost_thread-mt ) 
		#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops")
		#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops")
	ENDIF(WIN32)
ENDIF(CGAL_FOUND)

IF(WIN32)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /DCVT_MULTITHREAD ")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DCVT_MULTITHREAD ")
ELSE(WIN32)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DCVT_MULTITHREAD ")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCVT_MULTITHREAD ")
ENDIF(WIN32)



AUX_SOURCE_DIRECTORY(. SUB_SOURCES)
AUX_SOURCE_DIRECTORY(./algebra SOURCES)
AUX_SOURCE_DIRECTORY(./common SOURCES)
AUX_SOURCE_DIRECTORY(./combinatorics SOURCES)
AUX_SOURCE_DIRECTORY(./combinatorics/exact SOURCES)
SET(SOURCES ${SOURCES} ${SUB_SOURCES})
INCLUDE_DIRECTORIES(..)
ADD_EXECUTABLE(LpCVT ${SOURCES})

