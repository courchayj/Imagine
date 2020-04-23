// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// 3D affine apply
	template <typename T, int dim>
	class Affine {
	protected:
		typedef FMatrix<T,dim,dim> Matrix;
		typedef FVector<T,dim> Vector;
		Matrix _matrix; // Matrix
		Vector _vector; // Column vector

	public:
		// Constructors
		Affine() {}
		Affine(const Matrix& A, const Vector& b) : _matrix(A), _vector(b) {}
		Affine& operator = (const Affine& A) {
			_matrix = A._matrix;
			_vector = A._vector;
			return *this;
		}

		// Boolean comparison
		bool operator == (const Affine& B) const {
			return (_matrix == B._matrix && _vector == B._vector);
		}

		bool operator != (const Affine& B) const {
			return !(*this == B);
		}

		// Access to matrix and vector
		const Matrix& matrix() const { return _matrix; }
		Matrix& matrix() { return _matrix; }
		const Vector& vector() const { return _vector; }
		Vector& vector() { return _vector; }

		// Basic affine transformations
		static inline Affine Identity() {
			return Affine(Matrix::Identity(),Vector::Zero());
		}

		static inline Affine Translation(const Vector& v) {
			return Affine(Matrix::Identity(),v);
		}
		static inline Affine Translation(T tx,T ty) {
			return Affine(Matrix::Identity(),Vector(tx,ty));
		}
		static inline Affine Translation(T tx,T ty,T tz) {
			return Affine(Matrix::Identity(),Vector(tx,ty,tz));
		}

		static inline Affine Rotation(T theta){
			assert(dim==2);
			const T coth = cos(theta);
			const T sith = sin(theta);

			Matrix R;
			R(0,0) = coth;
			R(1,1) = coth;
			R(0,1) = sith;
			R(1,0) = -sith;

			return Affine(R,Vector::Zero());
		}

		static inline Affine Rotation(T alpha, T beta, T gamma) {
			assert(dim==3);
			const T coa = cos(alpha);
			const T sia = sin(alpha);
			const T cob = cos(beta);
			const T sib = sin(beta);
			const T cog = cos(gamma);
			const T sig = sin(gamma);

			Matrix _matrix;
			_matrix(0,0) = cog * cob;
			_matrix(0,1) = cog * sib * sia - sig * coa;
			_matrix(0,2) = cog * sib * coa + sig * sia;
			_matrix(1,0) = sig * cob;
			_matrix(1,1) = sig * sib * sia + cog * coa;
			_matrix(1,2) = sig * sib * coa - cog * sia;
			_matrix(2,0) = -sib;
			_matrix(2,1) = cob * sia;
			_matrix(2,2) = cob * coa;

			return Affine(_matrix,Vector::Zero());
		}


		static inline Affine Scaling(const Vector& v) {
			return Affine(Diagonal(v),Vector::Zero());
		}
		// nD isotropic
		static inline Affine Scaling(T s) { return Scaling(Vector(s)); }
		// 2D
		static inline Affine Scaling(T sx, T sy) { assert(dim==2); return Scaling(Vector(sx,sy)); }
		// 3D
		static inline Affine Scaling(T sx, T sy, T sz) { assert(dim==3); return Scaling(Vector(sx,sy,sz)); }

		// Composition
		Affine& operator *= (const Affine& B) {
			_vector += _matrix * B._vector;
			_matrix = _matrix * B._matrix;
			return *this;
		}

		Affine operator * (const Affine& B) const {
			Affine C;
			C._vector = _vector + _matrix * B._vector;
			C._matrix = _matrix * B._matrix;
			return C;
		}

		// Inverse affine apply
		friend Affine inverse(const Affine& A) {
			Affine B;
			B._matrix = inverse(A._matrix);
			B._vector = - B._matrix * A._vector;
			return B;
		}

		// Apply affine apply to a point
		Vector operator () (const Vector& v) const { return _matrix * v + _vector; }
		// (Faster?) specialized
		void apply(T x, T y, T& X, T& Y) const {
			assert(dim==2);
			X = _matrix(0,0) * x + _matrix(0,1) * y + _vector[0];
			Y = _matrix(1,0) * x + _matrix(1,1) * y + _vector[1];
		}
		// (Faster?) specialized
		void apply(T x, T y, T z, T& X, T& Y, T& Z) const {
			assert(dim==3);
			X = _matrix(0,0) * x + _matrix(0,1) * y + _matrix(0,2) * z + _vector[0];
			Y = _matrix(1,0) * x + _matrix(1,1) * y + _matrix(1,2) * z + _vector[1];
			Z = _matrix(2,0) * x + _matrix(2,1) * y + _matrix(2,2) * z + _vector[2];
		}

		// I/O
		friend void write (std::ostream& out, const Affine& A)  {
			write(out,A._matrix);
			write(out,A._vector);
		}

		friend void read (std::istream& in, Affine& A) {
			read(in,A._matrix);
			read(in,A._vector);
		}

		friend std::istream& operator >> (std::istream& in, Affine& A) {
			for (int i=0;i<dim;i++) {
				for (int j=0;j<dim;j++)
					in >> A._matrix(i,j);
				in >> A._vector[i];
			}
			return in;
		}

		friend std::ostream& operator << (std::ostream& out, const Affine& A) {
			for (int i=0;i<dim;i++) {
				for (int j=0;j<dim;j++)
					out << A._matrix(i,j) << " ";
				out << A._vector[i];
				out << std::endl;
			}
			return out;
		}
	};
	/// @}
}

