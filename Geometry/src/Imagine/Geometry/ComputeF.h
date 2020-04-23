// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{

	/// \brief Compute fundamental matrix.
	/// \details 7- and 8-point algorithms are implemented.
	template <typename T>
	class ComputeF {
	public:
		typedef FVector<T,2> Point2;
	protected:
		std::vector<Point2> ptsL, ptsR;
	public:
		ComputeF() {}
		void clear()
		{ ptsL.clear(); ptsR.clear(); }

		void addMatch(T xl, T yl, T xr, T yr)
		{ addMatch(Point2(xl,yl), Point2(xr,yr)); }
		void addMatch(const Point2& ptl, const Point2& ptr)
		{ ptsL.push_back(ptl); ptsR.push_back(ptr); }
		void addMatch(const FVector<T,4>& m)
		{ ptsL.push_back(Point2(m[0],m[1])); ptsR.push_back(Point2(m[2],m[3])); }

		/// Return number of solutions.
		template <class OutputIterator>
		int run(OutputIterator it)
		{
			if(ptsL.size() < 7)
				return 0;
			FMatrix<T,3,3> NL = normalize<T,2>(ptsL.begin(),ptsL.end());
			FMatrix<T,3,3> NR = normalize<T,2>(ptsR.begin(),ptsR.end());
			NL = transpose(NL);
			Matrix<T> A = fillMatrixEquations();
			Matrix<T> U,Vt; Vector<T> W;
			svd(A,U,W,Vt,true);
			int sols = 0;
			if(ptsL.size() >= 8) { // Force rank 2
				Vector<T> row=Vt.getRow(8);
				Matrix<T> f(row.data(),3,3); // Beware: shares memory. Don't kill row!
				svd(f,U,W,Vt);
				f.fill(T(0));
				f(0,0) = W[0]; f(1,1) = W[1];
				f = U*f*Vt;
				*it++ = Fundamental<T>(NL*FMatrix<T,3,3>(f.data())*NR);
				sols = 1;
			} else {
				FMatrix<T,3,3> f1(Vt.getRow(8).data()), f2(Vt.getRow(7).data());
				FMatrix<Polynomial<T,3>,3,3> m; // Matrix of polynomials
				for(int j=0; j < 3; j++)
					for(int i=0; i < 3; i++)
						m(i,j) = Polynomial<T,3>(f1(i,j)-f2(i,j),f2(i,j));
				Polynomial<T,3> P = det(m); // det(x F1 + (1-x)F2)
				T x[3];
				sols = P.cubicRoots(x[0],x[1],x[2]);
				for(int i=0; i < sols; i++)
					*it++ = Fundamental<T>(NL*(f1*x[i]+f2*(1-x[i]))*NR);
			}
			clear();
			return sols;
		}

	protected:
		Matrix<T> fillMatrixEquations() const
		{
			Matrix<T> A(int(ptsL.size()),9);
			for (int i=0; i < (int)ptsL.size(); i++) {
				A(i,0) = ptsL[i].x()*ptsR[i].x();
				A(i,1) = ptsL[i].y()*ptsR[i].x();
				A(i,2) =             ptsR[i].x();
				A(i,3) = ptsL[i].x()*ptsR[i].y();
				A(i,4) = ptsL[i].y()*ptsR[i].y();
				A(i,5) =             ptsR[i].y();
				A(i,6) = ptsL[i].x();
				A(i,7) = ptsL[i].y();
				A(i,8) = T(1);
			}
			return A;
		}
	};

	/// \brief Estimate fundamental matrix from a bunch of matches.
	/// \details This is just a wrapping around \c ComputeF to be used as functor
	/// for RANSAC.
	template<typename T>
	class FEstimator
	{
	public:
		template <typename DataIterator,int NumPoints,class OutputIterator>
		void operator()(const FArray<DataIterator,NumPoints>& pts,
			OutputIterator it) const
		{
			ComputeF<T> cmp;
			for(int i=0; i < pts.size(); i++)
				cmp.addMatch((*pts[i])[0], (*pts[i])[1], (*pts[i])[2], (*pts[i])[3]);
			cmp.run(it);
		}
	};

	/// \brief Functor evaluating distance of a point to its epipolar line.
	template <typename T>
	class FResidual
	{
	public:
		typedef enum {
			Left,     // Measure error in left image
			Right,    // Measure error in right image
			LeftRight // Measure error as max in both images
		} Direction;
	public:
		FResidual(Direction d=LeftRight)
			: dir(d) {}
		T operator()(const Fundamental<T>& F, const FVector<T,4>& match) const
		{
			typedef typename Fundamental<T>::Vector3 Vector3;
			float d1=0, d2=0;
			if(dir != Left) {
				Vector3 epi = F.rightLine(match[0],match[1]);
				d1 = epi*Vector3(match[2],match[3],T(1));
				if(d1 < 0) d1 = -d1;
				d1 /= std::sqrt(epi[0]*epi[0]+epi[1]*epi[1]);
			}
			if(dir != Right) {
				Vector3 epi = F.leftLine(match[2],match[3]);
				d2 = epi*Vector3(match[0],match[1],T(1));
				if(d2 < 0) d2 = -d2;
				d2 /= std::sqrt(epi[0]*epi[0]+epi[1]*epi[1]);
			}
			return std::max(d1,d2);
		}
	private:
		Direction dir;
	};

	/// @}
}

