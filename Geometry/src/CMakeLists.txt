# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(ImagineGeometry_MainHead Imagine/Geometry.h) 
SET(ImagineGeometry_Headers 
	Imagine/Geometry/Affine.h
	Imagine/Geometry/Camera.h
	Imagine/Geometry/Distortion.h
	Imagine/Geometry/Triangulation.h
	Imagine/Geometry/ComputeF.h
	Imagine/Geometry/TrifocalRansac.h
	Imagine/Geometry/Trifocal.h
	Imagine/Geometry/Polynomial.h
	Imagine/Geometry/Normalize.h
	Imagine/Geometry/Import.h
	Imagine/Geometry/Homography.h
	Imagine/Geometry/Graphics.h
	Imagine/Geometry/Fundamental.h
)

IF(IMAGINE_INSTALL)
INSTALL(FILES ${ImagineGeometry_MainHead} DESTINATION  include/Imagine)
INSTALL(FILES ${ImagineGeometry_Headers} DESTINATION  include/Imagine/Geometry)
ENDIF(IMAGINE_INSTALL)
