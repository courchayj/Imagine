// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

// SIFT. D. Lowe sources -> do not distribute!

namespace Imagine {
	/// \addtogroup Features
	/// @{
	
	/// Lowe SIFT descriptor
	typedef FeaturePoint<FVector<byte,128> > LoweSIFT;
	
	/// SIFT detector. Lowe's implementation.
	/// SIFT detector. Lowe's implementation.
	class LoweSIFTDetector:public FeatureDetector<LoweSIFT> {
		// First Octave Index. 
		int firstOctave;
		// Number of octaves.
		int numOctaves;	
		// Number of scales per octave. 
		int numScales;  
		// Max ratio of Hessian eigenvalues. 
		float edgeThresh;
		// Min contrast.
    	float peakThresh;
	public:
		/// Constructor.
		/// Constructor.
		/// 
		/// \dontinclude Features/test/test.cpp \skip main()
		/// \skipline Lowe sift, constructor
		LoweSIFTDetector() {
			firstOctave=-1;
			numOctaves=-1;
			numScales=3;
			edgeThresh=10.0f; peakThresh=0.04f;
		}
		/// Number of octaves.
		/// Sets number of octaves.  -1 = max. default=-1
		/// 
		/// \dontinclude Features/test/test.cpp \skip main()
		/// \skipline Lowe sift, number of octaves
		void setNumOctaves(int n) { numOctaves=n; }
		/// Number of scales per octave.
		/// Sets number of scales per octave. default=3
		/// 
		/// \dontinclude Features/test/test.cpp \skip main()
		/// \skipline Lowe sift, number of scales per octave
		void setNumScales(int n) { numScales=n; }
		/// First octave.
		/// Sets first octave. -1 = doubleImageSize. default=-1
		/// 
		/// \dontinclude Features/test/test.cpp \skip main()
		/// \skipline Lowe sift, first octave
		void setFirstOctave(int n) { firstOctave=n; }
		/// Max ratio of Hessian eigenvalues.
		/// Sets max ratio of Hessian eigenvalues. default=10
		/// 
		/// \dontinclude Features/test/test.cpp \skip main()
		/// \skipline Lowe sift, max ratio of Hessian eigenvalues
		void setEdgeThresh(float t) { edgeThresh=t; }
		/// Min contrast.
		/// Sets min contrast. Scale between 0 and 1. Default=0.04
		/// 
		/// \dontinclude Features/test/test.cpp \skip main()
		/// \skipline Lowe sift, min contrast
		void setPeakThresh(float t) { peakThresh=t; }
		
		// Implementation
		Array<LoweSIFT> run(const Image<byte>& I) const;
	};
	///@}

}
