# Imagine++ Libraries
# Copyright (C) Imagine
# For detailed information: http://imagine.enpc.fr/software

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(ImagineOptim_MainHead Imagine/Optim.h) 
SET(ImagineOptim_Headers 
	Imagine/Optim/block.h
	Imagine/Optim/energy.h
	Imagine/Optim/GCoptimization.h
	Imagine/Optim/GCoptimization_cpp.h
	Imagine/Optim/graph.h
	Imagine/Optim/graph_cpp.h
	Imagine/Optim/LinkedBlockList.h
	Imagine/Optim/LinkedBlockList_cpp.h
	Imagine/Optim/maxflow_cpp.h
	Imagine/Optim/Ransac.h
)

IF(IMAGINE_INSTALL)
INSTALL(FILES ${ImagineOptim_MainHead} DESTINATION  include/Imagine)
INSTALL(FILES ${ImagineOptim_Headers} DESTINATION  include/Imagine/Optim)
ENDIF(IMAGINE_INSTALL)
