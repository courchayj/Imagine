// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{

	template <typename T>
	class Homography : public FMatrix<T,3,3>
	{
	public:
		typedef FMatrix<T,3,3> Matrix3x3;
		typedef FVector<T,3> Vector3;
		typedef FVector<T,2> Vector2;

		Homography() { this->identity(); }
		Homography(const Matrix3x3& A) : Matrix3x3(A) {}
		Homography(const Camera<T>& cam1, const Camera<T>& cam2, const Vector3& point, const Vector3& normal)
		{			
			T nm = normal*point;
			Matrix3x3 p1n, p2n;

			for (int j=0; j<3; j++) for (int i=0; i<3; i++)
			{
				p1n(i,j) = cam1.vector()(i)*normal(j);
				p2n(i,j) = cam2.vector()(i)*normal(j);
			}

			Matrix3x3::operator = ( (cam2.matrix() + p2n / nm) * (cam1.matrix() + p1n / nm).inverse() );
		}

		void apply(T x, T y, T& X, T& Y) const
		{
			X = (*this)(0,0) * x + (*this)(0,1) * y + (*this)(0,2);
			Y = (*this)(1,0) * x + (*this)(1,1) * y + (*this)(1,2);
			T Z = (*this)(2,0) * x + (*this)(2,1) * y + (*this)(2,2);
			X /= Z;
			Y /= Z;
		}

		Vector2 apply(const Vector2& p) const
		{
			Vector2 P;
			apply(p.x(),p.y(),P.x(),P.y());
			return P;
		}
	};
	/// @}
}

