# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(ImagineFeaturesLib)

SET(HEADER Imagine/Features.h)
SET(SUBHEADERS
       Imagine/Features/Feat.h
       Imagine/Features/IO.h
       Imagine/Features/Match.h
       Imagine/Features/MSER.h
       Imagine/Features/SIFT.h
       Imagine/Features/LoweSIFT.h
)

ADD_LIBRARY(ImagineFeatures STATIC
       vl/generic.c
	    vl/host.c
	vl/imop.c 
	vl/sift.c
	vl/mser.c
	lowe/LoweKey.cpp
	lowe/LoweUtil.cpp
	MSER.cpp 
	SIFT.cpp
	LoweSIFT.cpp
       ${HEADER}
       ${SUBHEADERS}
       )
IF(WIN32)
	ImagineAppendProperty(ImagineFeatures DEBUG_POSTFIX "d")
ENDIF(WIN32)
ImagineUseModules(ImagineFeatures Images)

IF (IMAGINE_INSTALL)
	INSTALL(FILES ${SUBHEADERS} DESTINATION  include/Imagine/Features)
	INSTALL(FILES ${HEADER} DESTINATION  include/Imagine)
	INSTALL(TARGETS ImagineFeatures
            ARCHIVE DESTINATION  lib
            LIBRARY DESTINATION  lib
	)
	IF(WIN32)
		INSTALL(FILES 
			${ImagineFeaturesLib_BINARY_DIR}/Debug/ImagineFeaturesd.lib
			${ImagineFeaturesLib_BINARY_DIR}/Debug/ImagineFeaturesd.pdb
			${ImagineFeaturesLib_BINARY_DIR}/RelWithDebInfo/ImagineFeatures.pdb
			DESTINATION  lib)
	ENDIF(WIN32)
ENDIF (IMAGINE_INSTALL)
