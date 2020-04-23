# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

# Loading 
IF(NOT ImagineOptim_FOUND)
SET(ImagineOptim_FOUND true)

IF (IMAGINE_FROM_SRC)
	SET(d  $ENV{IMAGINEPP_SRC}/Imagine-labo/Optim/src)
		INCLUDE_DIRECTORIES(BEFORE ${d})
ELSE (IMAGINE_FROM_SRC)
	SET(d  $ENV{IMAGINEPP_ROOT}/include)
		INCLUDE_DIRECTORIES(AFTER ${d})
ENDIF (IMAGINE_FROM_SRC)
SET(IMAGINE_OPTIM_HEADERS
	${d}/Imagine/Optim.h 
	${d}/Imagine/Optim/GCoptimization.h
	${d}/Imagine/Optim/Graph.h
	${d}/Imagine/Optim/Ransac.h
	)
FILE(TO_CMAKE_PATH "${IMAGINE_OPTIM_HEADERS}" IMAGINE_OPTIM_HEADERS)
ENDIF(NOT ImagineOptim_FOUND)

# Using Optim
IF (NOT ${Optim_Proj} STREQUAL "_PRELOAD_")
	ImagineUseModules(${Optim_Proj} Common)
ENDIF (NOT ${Optim_Proj} STREQUAL "_PRELOAD_")
