# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

# Loading 
IF(NOT ImagineGeometry_FOUND)
SET(ImagineGeometry_FOUND true)

IF (IMAGINE_FROM_SRC)
	SET(d  $ENV{IMAGINEPP_SRC}/Imagine-labo/Geometry/src)
		INCLUDE_DIRECTORIES(BEFORE ${d})
ELSE (IMAGINE_FROM_SRC)
	SET(d  $ENV{IMAGINEPP_ROOT}/include)
		INCLUDE_DIRECTORIES(AFTER ${d})
ENDIF (IMAGINE_FROM_SRC)
SET(IMAGINE_GEOMETRY_HEADERS
	${d}/Imagine/Geometry.h 
	${d}/Imagine/Geometry/Affine.h
	${d}/Imagine/Geometry/Camera.h
	${d}/Imagine/Geometry/Distortion.h
	${d}/Imagine/Geometry/Triangulation.h
	${d}/Imagine/Geometry/ComputeF.h
	${d}/Imagine/Geometry/TrifocalRansac.h
	${d}/Imagine/Geometry/Trifocal.h
	${d}/Imagine/Geometry/Polynomial.h
	${d}/Imagine/Geometry/Normalize.h
	${d}/Imagine/Geometry/Import.h
	${d}/Imagine/Geometry/Homography.h
	${d}/Imagine/Geometry/Graphics.h
	${d}/Imagine/Geometry/Fundamental.h
)
FILE(TO_CMAKE_PATH "${IMAGINE_GEOMETRY_HEADERS}" IMAGINE_GEOMETRY_HEADERS)

ENDIF(NOT ImagineGeometry_FOUND)

# Using Geometry
IF (NOT ${Geometry_Proj} STREQUAL "_PRELOAD_")
	ImagineUseModules(${Geometry_Proj} Common)
	ImagineUseModules(${Geometry_Proj} Graphics)
	ImagineUseModules(${Geometry_Proj} Images)
	ImagineUseModules(${Geometry_Proj} LinAlg)
	ImagineUseModules(${Geometry_Proj} Optim)
ENDIF (NOT ${Geometry_Proj} STREQUAL "_PRELOAD_")
