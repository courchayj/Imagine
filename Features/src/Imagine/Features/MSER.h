// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

//MSER + SIFT descriptors. From VLFeat

namespace Imagine {
	/// \addtogroup Features
	/// @{

	/// MSER descriptor
	typedef FeaturePoint<FVector<byte,128> > MSER;
	
	/// MSER detector. VLFeat implementation.
	/// MSER detector. VLFeat implementation.
	class MSERDetector:public FeatureDetector<MSER> {
	public:
		/// Constructor.
		/// Constructor.
		/// 
		/// \dontinclude Features/test/test.cpp \skip main()
		/// \skipline MSER, constructor
		MSERDetector() {		
		}
		// Implementation
		Array<MSER> run(const Image<byte>& I) const;
	};
	///@}

}
