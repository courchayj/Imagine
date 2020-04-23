// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{


	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Class for a perspective camera
	// The normalization of the camera matrix is enforced: the norm of the third row is equal to 1
	template <typename T>
	class Camera : protected Affine<T,3> {

	public:
		typedef Affine<T,3> Affine3;
		typedef FMatrix<T,3,3> Matrix3x3;
		typedef FVector<T,3> Vector3;
		typedef FVector<T,2> Vector2;

	protected:
		// Parametres extrinseques
		Affine3 _extrinsic;
		Vector3 _center;

		// Parametres intrinseques
		Vector2 _principal_point, _pixel_size;
		T _skew;

		// Pour la back projection
		Affine3 _inverse;

		// Radial distortion parameters
		Vector2 _distortion_center;
		Vector3 _distortion_coeffs;
		bool _distortion_direct;

	public:
		// Constructeur par defaut
		Camera(): Affine3(Affine3::Identity()) { updateParameters(); resetDistortion(); }

		// Constructor from a 3x4 matrix
		Camera(const Affine3& A) : Affine3(A) { updateParameters(); resetDistortion(); }
		Camera(const Matrix3x3& A, const Vector3& b) : Affine3(A,b) { updateParameters(); resetDistortion(); }

		// Constructor from extrinsic and intrinsic parameters
		// Caution : the extrinsic parameters must be a rigid transformation,
		// Otherwise, the camera representation may be inconsistent
		Camera(const Affine3& extrinsic, const Vector2& principalPoint, const Vector2& pixelSize, T skew=T(0))
			: _extrinsic(extrinsic), _principal_point(principalPoint), _pixel_size(pixelSize),_skew(skew) {
				updateMatrix();
				resetDistortion();
		}

		Camera(const Affine3& extrinsic, T principalX, T principalY, T pixelWidth, T pixelHeight, T skew=T(0))
			: _extrinsic(extrinsic), _principal_point(principalX,principalY), _pixel_size(pixelWidth,pixelHeight) ,_skew(skew){
				updateMatrix();
				resetDistortion();
		}

		// Constructeur par recopie
		Camera& operator = (const Camera& cam) {
			Affine3::operator = (cam);
			_extrinsic = cam._extrinsic;
			_principal_point = cam._principal_point;
			_pixel_size = cam._pixel_size;
			_skew=cam._skew;
			_center = cam._center;
			_distortion_center = cam._distortion_center;
			_distortion_coeffs = cam._distortion_coeffs;
			_distortion_direct = cam._distortion_direct;
			updateCenter();
			return *this;
		}

		Camera& operator = (const Affine3& A) {
			Affine3::operator = (A);
			updateParameters();
			resetDistortion();
			return *this;
		}

		// Boolean comparison
		bool operator == (const Camera& cam) const {
			return (
				Affine3::operator == (cam.affine()) &&
				_distortion_center == cam.distortionCenter() &&
				_distortion_coeffs == cam.distortionCoeffs() &&
				_distortion_direct == cam.distortionDirect()
				);
		}

		bool operator != (const Camera& cam) const {
			return !(*this == cam);
		}

		// Access to matrix
		const Affine3& affine() const { return *this; }
		const Matrix3x3& matrix() const { return Affine3::matrix(); }
		const Vector3& vector() const { return Affine3::vector(); }

		// Access to extrinsic parameters
		const Affine3& extrinsic() const { return _extrinsic; }

		// Caution : the extrinsic parameters must be a rigid transformation,
		// Otherwise, the camera representation may be inconsistent
		void setExtrinsic(const Affine3& A) {
			_extrinsic = A;
			updateMatrix();
		}


		// Access to intrinsic parameters
		const Vector2& principalPoint() const { return _principal_point; }
		T principalX() const { return _principal_point.x(); }
		T principalY() const { return _principal_point.y(); }
		const Vector2& pixelSize() const { return _pixel_size; }
		T pixelWidth() const { return _pixel_size.x(); }
		T pixelHeight() const { return _pixel_size.y(); }
		T skew() const { return _skew; }


		Camera& setIntrinsic(T principalX, T principalY, T pixelWidth, T pixelHeight, T skew=T(0)) {
			return setIntrinsic(Vector2(principalX,principalY), Vector2(pixelWidth,pixelHeight),skew);
		}

		Camera& setIntrinsic(const Vector2& principalPoint, const Vector2& pixelSize, T skew=T(0)) {
			_principal_point = principalPoint;
			_pixel_size = pixelSize;
			_skew=skew;
			return updateMatrix();
		}


		// Access to distortion parameters
		void resetDistortion()
		{
			_distortion_center = _principal_point;
			_distortion_coeffs.fill(0.);
			_distortion_direct = true;
		}

		bool hasDistortion() const 
		{
			return (_distortion_coeffs != Vector3(T(0)));
		}

		const Vector2& distortionCenter() const { return _distortion_center; }
		Vector2& distortionCenter() { return _distortion_center; }
		const Vector3& distortionCoeffs() const { return _distortion_coeffs; }
		Vector3& distortionCoeffs() { return _distortion_coeffs; }
		bool distortionDirect() const { return _distortion_direct; }
		bool& distortionDirect() { return _distortion_direct; }


		// Centre de la camera
		const Vector3& center() const { return _center; }		

		// Direction de la camera
		Vector3 direction() const { return this->_matrix.getRow(2); }

		// Profondeur d'un point
		T depth(T x, T y, T z) const
		{
			return this->matrix()(2,0) * x + this->matrix()(2,1) * y + this->matrix()(2,2) * z + this->vector()[2];
		}

		// Projection d'un point
		inline void projection(T x, T y, T z, T& X, T& Y) const {
			T u,v,w;
			this->apply(x,y,z,u,v,w);
			X = u/(w+1e-12);
			Y = v/(w+1e-12);
		}
		// Projection d'un point
		inline FVector<T,2> projection(FVector<T,3> M) const {
			T X,Y;
			this->projection(M[0],M[1],M[2],X,Y);
			return FVector<T,2>(X,Y);
		}

		// Projection d'un point et profondeur associee
		inline void projection(T x, T y, T z, T& X, T& Y, T& Z) const {
			T u,v;
			this->apply(x,y,z,u,v,Z);
			X = u/Z;
			Y = v/Z;
		}

		void backProjection(T X, T Y, T Z, T& x, T& y, T& z) const {
			_inverse.apply(X*Z,Y*Z,Z,x,y,z);
		}


		// Add or remove radial distortion to pixel positions
		void addDistortion(T xu, T yu, T& xd, T& yd) const
		{
			xu = ( xu - distortionCenter().x() ) / pixelWidth();
			yu = ( yu - distortionCenter().y() ) / pixelHeight();
			T r = std::sqrt(xu*xu+yu*yu);
			T lambda = (r==0) ? T(1) : (distortionDirect() ? radialDirect(r) : radialInverse(r)) / r;
			xd = distortionCenter().x() + lambda * xu * pixelWidth();
			yd = distortionCenter().y() + lambda * yu * pixelHeight();
		}

		void removeDistortion(T xd, T yd, T& xu, T& yu) const
		{
			xd = ( xd - distortionCenter().x() ) / pixelWidth();
			yd = ( yd - distortionCenter().x() ) / pixelHeight();
			T r = sqrt(xd*xd+yd*yd);
			T lambda = (r==0) ? T(1) : (distortionDirect() ? radialInverse(r) : radialDirect(r)) / r;
			xu = distortionCenter().x() + lambda * xd * pixelWidth();
			yu = distortionCenter().y() + lambda * yd * pixelHeight();
		}


		// Composition by an affine transfo
		Camera operator * (const Affine3& B) const {
			Camera C(*this);
			return (C *= B);
		}

		Camera& operator *= (const Affine3& B) {
			Affine3::operator *= (B);
			return updateParameters();
		}

		// Application d'un facteur d'echelle a l'espace ou a l'image
		Camera& scaleWorld(T f) {
			return scaleWorld(f,f,f);
		}
		Camera& scaleWorld(T fx, T fy, T fz) {
			this->_matrix.setCol(0, this->_matrix.getCol(0) * fx);
			this->_matrix.setCol(1, this->_matrix.getCol(1) * fy);
			this->_matrix.setCol(2, this->_matrix.getCol(2) * fz);
			return updateParameters();
		}
		Camera& scaleImage(T f) {
			return scaleImage(f,f);
		}
		Camera& scaleImage(T fx, T fy) {
			this->_matrix.setRow(0, this->_matrix.getRow(0) * fx);
			this->_matrix.setRow(1, this->_matrix.getRow(1) * fy);
			this->_vector[0] *= fx;
			this->_vector[1] *= fy;
			_distortion_center.x() *= fx;
			_distortion_center.y() *= fy;
			return updateParameters();
		}

		// Application d'une translation a l'image
		Camera& translateImage(T tx, T ty)
		{
			this->_matrix.setRow(0, this->_matrix.getRow(0) + tx * this->_matrix.getRow(2));
			this->_matrix.setRow(1, this->_matrix.getRow(1) + ty * this->_matrix.getRow(2));
			this->_vector[0] += tx * this->_vector[2];
			this->_vector[1] += ty * this->_vector[2];
			_distortion_center += Vector2(tx,ty);
			return updateParameters();
		}

		// Application d'une translation au monde
		Camera& translateWorld(T tx, T ty, T tz)
		{
			this->_vector += this->_matrix * Vector3(tx,ty,tz);
			return updateParameters();
		}

		// Epipole d'une autre camera dans cette camera (en coordonees projective car il peut etre a l'infini)
		Vector3 epipole(const Camera& cam) const {
			return (*this)(cam.center());
		}

		// Matrice fondamentale d'une paire de cameras
		Matrix3x3 fundamental(const Camera& cam) const {
			return Matrix3x3::CrossProd(epipole(cam)) * this->_matrix * inverse(cam._matrix);
		}

		// Entrees sorties
		friend std::istream& operator >> (std::istream& in, Camera& cam) {
			in >> (Affine3&) cam;
			cam.updateParameters();
			cam.resetDistortion();
			return in;
		}
		friend std::ostream& operator << (std::ostream& o, const Camera& cam) {
			o << (const Affine3&) cam;
			return o;
		}
		friend void write (std::ostream& o, const Camera& cam)  { 
			write(o,(const Affine3&)cam); 
		}
		friend void read(std::istream& in, Camera& cam) {
			read(in,(Affine3&)cam); cam.updateParameters(); cam.resetDistortion(); 
		}



		// Normalisation de la matrice de camera : la troisieme ligne doit etre de norme 1
		void normalize() {
			const T n = norm(this->_matrix.getRow(2));
			this->_matrix /= n;
			this->_vector /= n;
		}
		void givens(const Matrix3x3& A,Matrix3x3& R,Matrix3x3& Q) { // 3x3 RQ decomp + diag R > 0
			// Beware: R: intrinsic params (cam.K); Q: rotation (cam.R)
			T cs, sn, a, b;
			Matrix3x3 Qx,Qy,Qz;
			R=A;
			a=R(2,1);b=R(2,2);
			cs = - b / std::sqrt(a * a + b * b);
			sn = a / std::sqrt(a * a + b * b);
			Qx=Matrix3x3::Identity();
			Qx(1,1)=Qx(2,2)=cs;
			Qx(1,2)=-sn;Qx(2,1)=sn;
			R=R*Qx;
			a=R(2,0);b=R(2,2);
			cs = b / std::sqrt(a * a + b * b);
			sn = a / std::sqrt(a * a + b * b);
			Qy=Matrix3x3::Identity();
			Qy(0,0)=Qy(2,2)=cs;
			Qy(0,2)=sn;Qy(2,0)=-sn;
			R=R*Qy;
			a=R(1,0);b=R(1,1);
			cs = -b / std::sqrt(a * a + b * b);
			sn = a / std::sqrt(a * a + b * b);
			Qz=Matrix3x3::Identity();
			Qz(0,0)=Qz(1,1)=cs;
			Qz(0,1)=-sn;Qz(1,0)=sn;
			R=R*Qz;
			Q=transpose(Qx*Qy*Qz);

			//Rectify matrix if needed
			if (R(0,0)<0) {
				R(0,0)=-R(0,0);
				Q.setRow(0,-Q.getRow(0));
			}
			if (R(1,1)<0) {
				R(0,1)=-R(0,1);
				R(1,1)=-R(1,1);
				Q.setRow(1,-Q.getRow(1));
			}
		}


		// Mise a jour des parametres a partir de la matrice

	protected:
		Camera& updateParameters() {
			normalize();
			Matrix3x3 K,R;
			givens(this->_matrix,K,R);
			_extrinsic=Affine3(R,inverse(K)*this->_vector);
			_pixel_size = Vector2(K(0,0),K(1,1));
			_skew=K(0,1);
			_principal_point = Vector2(K(0,2),K(1,2));
			return updateCenter();
		}

		// Mise a jour de la matrice a partir des parametres
		Camera& updateMatrix() {
			Matrix3x3 K;
			K(0,0)=_pixel_size[0];K(0,1)=_skew;K(0,2)=_principal_point[0];
			K(1,0)=T(0);K(1,1)=_pixel_size[1];K(1,2)=_principal_point[1];
			K(2,0)=K(2,1)=T(0);K(2,2)=T(1);
			Affine3::operator=(Affine3(K*_extrinsic.matrix(),K*_extrinsic.vector()));
			return updateCenter();
		}

		Camera& updateCenter() {
			_inverse = inverse(*this);
			_center = - _inverse.matrix() * this->vector();
			return *this;
		}

		T radialDirect(T r) const
		{
			T r2 = r*r;
			return r * (T(1) + r2 * (_distortion_coeffs[0] + r2 * (_distortion_coeffs[1] + r2 * _distortion_coeffs[2])));
		}

		T radialInverse(T r) const
		{
			// We use dichotomy to inverse the radial function
			T r1=r*T(0.9), r2=r*T(1.1); // Enough?
			bool s1=(r<radialDirect(r1));
			bool s2=(r<radialDirect(r2));
			if (s1==s2)
			{
				std::cerr << "failed to inverse radial dirtortion" << std::endl;
				return r;
			}
			for (int k=0;k<10;k++)
			{
				T rmid=(r1+r2)/T(2);
				bool smid=(r<radialDirect(rmid));
				if (s1==smid)
					r1=rmid;
				else
					r2=rmid;
			}
			return r1;
		}		
	};
	/// @}
}

