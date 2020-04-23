# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF (IMAGINE_INSTALL)
	INSTALL(FILES CMake/UseImagineGeometry.cmake DESTINATION  CMake)
ENDIF (IMAGINE_INSTALL)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
