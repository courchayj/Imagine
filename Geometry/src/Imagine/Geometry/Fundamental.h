// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{


	/// Class for a fundamental matrix.
	template <typename T >
	class Fundamental  {

	public:		
		typedef FMatrix<T,3,3> Matrix3x3;
		typedef FVector<T,3> Vector3;

	public:
		Fundamental() {}
		Fundamental(const Matrix3x3& F): _matrix(F) {}

		// Constructeur par recopie
		Fundamental& operator = (const Fundamental& fond) 
		{			
			_matrix = fond._matrix;
			return *this;
		}

		// Access to matrix, camera, epipole
		const Matrix3x3& matrix() const { return _matrix; }
		Vector3 leftLine(const Vector3& ptRight) const
		{ return _matrix*ptRight; }
		Vector3 rightLine(const Vector3& ptLeft) const
		{ return tmult(_matrix,ptLeft); }

		Vector3 leftLine(T x, T y) const
		{ return leftLine( Vector3(x,y,T(1)) ); }
		Vector3 rightLine(T x, T y) const
		{ return rightLine( Vector3(x,y,T(1)) ); }

		/// Epipole in left image.
		Vector3 leftEpipole()  const
		{
			Matrix<T> F(3,3), U,Vt;
			F(0,0) = _matrix(0,0);
			F(0,1) = _matrix(0,1);
			F(0,2) = _matrix(0,2);
			F(1,0) = _matrix(1,0);
			F(1,1) = _matrix(1,1);
			F(1,2) = _matrix(1,2);
			F(2,0) = _matrix(2,0);
			F(2,1) = _matrix(2,1);
			F(2,2) = _matrix(2,2);
			Vector<T> W;
			F.svda(U,W,Vt);
			return Vector3(U.getCol(2).data());
		}
		/// Epipole in right image.
		Vector3 rightEpipole() const
		{
			Matrix<T> F(3,3), U,Vt;
			F(0,0) = _matrix(0,0);
			F(0,1) = _matrix(0,1);
			F(0,2) = _matrix(0,2);
			F(1,0) = _matrix(1,0);
			F(1,1) = _matrix(1,1);
			F(1,2) = _matrix(1,2);
			F(2,0) = _matrix(2,0);
			F(2,1) = _matrix(2,1);
			F(2,2) = _matrix(2,2);
			Vector<T> W;
			F.svda(U,W,Vt);
			return Vector3(Vt.getRow(2).data());
		}
	protected:
		Matrix3x3 _matrix;
	};
	/// @}
}

