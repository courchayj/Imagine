// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine {
	/// \addtogroup Features
	/// @{

	/// Two best matches.
	/// Find the two best matches for a given FP. Note: distance to second best match is available. Not its index.
	/// \tparam T FP type
	/// \param F1,F2 FP arrays to match
	/// \param i1 index of FP of F1 to match
	/// \param bi2 index of best match in F2
	/// \param bd distance between F1[i1].desc and F2[bi2].desc
	/// \param sbd distance between F1[i1].desc and descriptor of second best match
	/// \param F pointer to fundamental matrix (optional)
	/// \param epiDist max distance to epipolar of F1[i1].pos
	/// 
	/// \dontinclude Features/test/test.cpp \skip matching()
	/// \skipline two best matches
	/// \until ...
	template <typename T>
	void find2Best(const Array<T>& F1,const Array<T>& F2,
		size_t i1,size_t& bi2,double& bd,double& sbd,FMatrix<double,3,3>* F=0,double epiDist=0
		) 
	{
		FVector<double,3> m1(F1[i1].pos[0],F1[i1].pos[1],1);
		FVector<double,3> l;
		if (F)
			l=(*F)*m1;

		bi2=-1;
		bd=sbd=std::numeric_limits<double>::max();
		for (size_t i2=0;i2<F2.size();i2++) {
			if (F) {
				FVector<double,3> m2(F2[i2].pos[0],F2[i2].pos[1],1);
				if (abs(l*m2)/std::sqrt(l[0]*l[0]+l[1]*l[1])>epiDist)
					continue;
			}

			double d=dist(F1[i1].desc,F2[i2].desc);
			if (d<bd) {
				sbd=bd;
				bd=d;
				bi2=i2;
			} else if (d<sbd) {
				sbd=d;			
			}
		}
		
	}
	
	/// Alias to list of pairs of indices
	typedef std::list<std::pair<size_t,size_t> > MatchList;

	/// Lowe's match.
	/// Matching based on Lowe's recommendations: ratio between best and second best matches. Naive impl. without ANN
	/// \tparam T FP type
	/// \param F1,F2 FP arrays to match
	/// \param ratio First/Second ratio
	/// \param sym symmetric test: match OK if OK both ways
	/// \param descDist maximum descriptor distance (0=none)
	/// \param F pointer to fundamental matrix (optional)
	/// \param epiDist max distance to epipolar of F1[i1].pos
	/// \return list of matches
	/// 
	/// \dontinclude Features/test/test.cpp \skip matching()
	/// \skipline Lowe's match
	/// \until ...
	template <typename T>
	MatchList loweMatch(
		const Array<T>& F1,
		const Array<T>& F2,
		double ratio=0.6,
		bool sym=false,
		double descDist=0,	
		FMatrix<double,3,3>* F=0,double epiDist=0
		)
	{	
			MatchList L;
			FMatrix<double,3,3> FTrans; // F transpose
			FMatrix<double,3,3>* FT=0;	// pointer (or 0)
			if (F) {
				FTrans=transpose(*F);
				FT=&FTrans;
			}
			for (size_t i1=0;i1<F1.size();i1++) {
				size_t bi2;
				double bd,sbd;
				find2Best(F1,F2,i1,bi2,bd,sbd,F,epiDist);
				if (bi2!=-1 && bd/sbd < ratio && (descDist==0 || bd <descDist)) {
					if (!sym) 
						L.push_back(std::pair<size_t,size_t>(i1,bi2));
					else {
						size_t bi1;
						find2Best(F2,F1,bi2,bi1,bd,sbd,FT,epiDist);
						if (bd/sbd<ratio && bi1==i1)	// i1 -> bi2 -> bi1 with bi1==i1
							L.push_back(std::pair<size_t,size_t>(i1,bi2));
					}
				}
			}
			return L;
	}
	///@}

}
