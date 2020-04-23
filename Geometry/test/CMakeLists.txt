# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

PROJECT(ImagineGeometryTests)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{IMAGINEPP_ROOT}/CMake)
FIND_PACKAGE(Imagine)

ImaginePreLoadModules(Geometry)
ImagineAddExecutable(ImagineGeometryExample example.cpp ${IMAGINE_GEOMETRY_HEADERS})
ImagineUseModules(ImagineGeometryExample Geometry)

IF(IMAGINE_INSTALL)
	INSTALL(FILES CMakeLists.txt example.cpp DESTINATION test/Geometry)
ENDIF(IMAGINE_INSTALL)

