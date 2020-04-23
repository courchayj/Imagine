# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

PROJECT(ImagineOptimTests)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{IMAGINEPP_ROOT}/CMake)
FIND_PACKAGE(Imagine)

ImaginePreLoadModules(Optim)
ImagineAddExecutable(ImagineOptimExample example.cpp ${IMAGINE_OPTIM_HEADERS})
ImagineUseModules(ImagineOptimExample Optim)

IF(IMAGINE_INSTALL)
	INSTALL(FILES CMakeLists.txt example.cpp DESTINATION test/Optim)
ENDIF(IMAGINE_INSTALL)


