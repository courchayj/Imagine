# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ImagineAddSubDirectory(Features)
ImagineAddSubDirectory(Optim)
ImagineAddSubDirectory(Geometry)
ImagineAddSubDirectory(LevelSets)
IF(WIN32)
	ImagineAddSubDirectory(GPU)
ENDIF(WIN32)
