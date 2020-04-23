// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{

	/// \brief A polynomial in one variable of degree <= \a MaxDegree.
	template <class T, int MaxDegree>
	class Polynomial : public FArray<T, MaxDegree+1>
	{
	public:
		Polynomial(): FArray<T,MaxDegree+1>()
		{ FArray<T,MaxDegree+1>::fill(0); }
		Polynomial(T a): FArray<T,MaxDegree+1>()
		{ FArray<T,MaxDegree+1>::fill(0); (*this)[0] = a; }
		Polynomial(T aX, T b): FArray<T,MaxDegree+1>()
		{ FArray<T,MaxDegree+1>::fill(0); (*this)[0] = b; (*this)[1] = aX; }

		T evaluate(T x) const
		{
			T y=T(0);
			for(int i=MaxDegree; i >= 0; i--)
				y = (*this)[i] + x*y;
			return y;
		}

		Polynomial<T,MaxDegree> operator+(const Polynomial<T,MaxDegree>& P) const
		{
			Polynomial<T,MaxDegree> Q=*this;
			Q += P;
			return Q;
		}
		Polynomial<T,MaxDegree> operator-(const Polynomial<T,MaxDegree>& P) const
		{
			Polynomial<T,MaxDegree> Q=*this;
			Q -= P;
			return Q;
		}
		void operator+=(const Polynomial<T,MaxDegree>& P)
		{ for(int i=0; i <= MaxDegree; i++) (*this)[i] += P[i]; }
		void operator-=(const Polynomial<T,MaxDegree>& P)
		{ for(int i=0; i <= MaxDegree; i++) (*this)[i] -= P[i]; }

		Polynomial<T,MaxDegree> operator*(const Polynomial<T,MaxDegree>& P) const
		{
			Polynomial<T,MaxDegree> Q;
			for(int i=0; i <= MaxDegree; i++)
				for(int j=0; i+j <= MaxDegree; j++)
					Q[i+j] += (*this)[i]*P[j];
			return Q;
		}

		int affineRoots(T& x0) const
		{
			if((*this)[1] == T(0))
				return 0;
			x0 = -(*this)[0] / (*this)[1];
			return 1;
		}

		int quadraticRoots(T& x0, T& x1) const
		{
			if((*this)[2] == T(0))
				return affineRoots(x0);
			T d = (*this)[1] - 4*(*this)[0]*(*this)[2];
			if(d < T(0))
				return 0;
			d = std::sqrt(d);
			if(d == T(0)) {
				x0 = -(*this)[1]/(2*(*this)[2]);
				return 1;
			}
			x0 = ((*this)[1]>=T(0))? -((*this)[1]+d): -(*this)[1]+d;
			x1 = 2*(*this)[0]/x0;
			x0 /= (2*(*this)[2]);
			return 2;
		}

		int cubicRoots(T& x0, T& x1, T& x2) const
		{
			if((*this)[3] == T(0))
				return quadraticRoots(x0, x1);
			T a1 = (*this)[2] / (*this)[3];
			T a2 = (*this)[1] / (*this)[3];
			T a3 = (*this)[0] / (*this)[3];

			T Q = (a1 * a1 - 3 * a2) / 9;
			T R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
			T Q3 = Q * Q * Q;
			T d = Q3 - R * R;

			if (d >= 0) { // Three real roots
				T theta = std::acos(R / std::sqrt(Q3));
				T sqrtQ = std::sqrt(Q);
				x0 = -2 * sqrtQ * std::cos( theta             / 3) - a1 / 3;
				x1 = -2 * sqrtQ * std::cos((theta + 2 * T(M_PI)) / 3) - a1 / 3;
				x2 = -2 * sqrtQ * std::cos((theta + 4 * T(M_PI)) / 3) - a1 / 3;
				return 3;
			}
			T e = std::pow(std::sqrt(-d) + std::fabs(R), T(1. / 3.));
			if (R > 0)
				e = -e;
			x0 = (e + Q / e) - a1 / 3;
			return 1;
		}
	};

	/// @}
}

