# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

# Loading ImagineFeatures
IF(NOT ImagineFeatures_FOUND)
SET(ImagineFeatures_FOUND true)

IF (IMAGINE_FROM_SRC)
	SET(d $ENV{IMAGINEPP_SRC}/Imagine-labo/Features/src)
		INCLUDE_DIRECTORIES(BEFORE ${d})
ELSE (IMAGINE_FROM_SRC)
	SET(d  $ENV{IMAGINEPP_ROOT}/include)
		INCLUDE_DIRECTORIES(AFTER ${d})
ENDIF (IMAGINE_FROM_SRC)
SET(IMAGINE_FEATURES_HEADERS
	${d}/Imagine/Features.h
	${d}/Imagine/Features/Feat.h
	${d}/Imagine/Features/IO.h
	${d}/Imagine/Features/Match.h
	${d}/Imagine/Features/MSER.h
	${d}/Imagine/Features/SIFT.h
	${d}/Imagine/Features/LoweSIFT.h
)
FILE(TO_CMAKE_PATH "${IMAGINE_FEATURES_HEADERS}" IMAGINE_FEATURES_HEADERS)
ENDIF(NOT ImagineFeatures_FOUND)

# Using Features 
IF (NOT ${Features_Proj} STREQUAL "_PRELOAD_")
	ImagineUseModules(${Features_Proj} Images)
	IF (IMAGINE_FROM_SRC)
		# Refers to project name
		TARGET_LINK_LIBRARIES(${Features_Proj} ImagineFeatures)
	ELSE (IMAGINE_FROM_SRC)
		# Refers to .lib names 
		FILE(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}" d)
		ImagineAddLinkPath(${Features_Proj} "${d}/lib")
		IF(WIN32)
			TARGET_LINK_LIBRARIES(${Features_Proj} debug ImagineFeaturesd)
			TARGET_LINK_LIBRARIES(${Features_Proj} optimized ImagineFeatures)
		ELSE(WIN32)
			TARGET_LINK_LIBRARIES(${Features_Proj} ImagineFeatures)
		ENDIF(WIN32)
	ENDIF (IMAGINE_FROM_SRC)
ENDIF (NOT ${Features_Proj} STREQUAL "_PRELOAD_")
