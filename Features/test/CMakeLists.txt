# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

PROJECT(ImagineFeaturesTests)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{IMAGINEPP_ROOT}/CMake)
FIND_PACKAGE(Imagine)

ImaginePreLoadModules(Features)
ImagineAddExecutable(ImagineFeaturesTest test.cpp ${IMAGINE_FEATURES_HEADERS})
ImagineUseModules(ImagineFeaturesTest Features)

ImagineAddExecutable(ImagineFeaturesExample example.cpp)
ImagineUseModules(ImagineFeaturesExample Features)

IF(IMAGINE_INSTALL)
	INSTALL(FILES CMakeLists.txt test.cpp example.cpp toys.0000.png toys.0200.png toys_F_0000_0200.txt DESTINATION test/Features)
ENDIF(IMAGINE_INSTALL)

