// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{


	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Class for a fondamental matrix

	template <typename T , int nbPoints,typename DataInputIterator>
	class Trifocal  {

	public:
		typedef FMatrix<T,3,3> Matrix3x3;
		typedef FVector<T,3> Vector3;
		typedef Camera<T> Cam;
		typedef Vector<T> Vectr;
		typedef Matrix<T> Matr;


	protected:

		Cam _camera1;
		Cam _camera2;
		Cam _camera3;
	public:
		// Constructeur par defaut
		Trifocal() {}

		//Generate Trifocal Tenseur from 7 pts 3 images at least
		Trifocal(  const FArray<DataInputIterator,nbPoints>& points)
		{
			Matr T1(3,3),T2(3,3),T3(3,3);
			Matr X1_pts(3,nbPoints), X2_pts(3,nbPoints),X3_pts(3,nbPoints);
			X1_pts.fill(T(1));
			X2_pts.fill(T(1));
			X3_pts.fill(T(1));
			//convert points in matrix
			for (int i = 0; i<nbPoints; i++)
			{				
				X1_pts(0,i) = (*points[i])[0];
				X1_pts(1,i) = (*points[i])[1];
				X2_pts(0,i) = (*points[i])[2];
				X2_pts(1,i) = (*points[i])[3];
				X3_pts(0,i) = (*points[i])[4];
				X3_pts(1,i) = (*points[i])[5];
			}

			T1 = normalize2D(X1_pts);
			T2 = normalize2D(X2_pts);
			T3 = normalize2D(X3_pts);
			X1_pts = T1 * X1_pts;
			X2_pts = T2 * X2_pts;
			X3_pts = T3 * X3_pts;

			computeMatrix(T1,T2,T3,X1_pts,X2_pts,X3_pts);

		}

		//Generate Trifocal Tenseur T213 from 2 fundametal matrices F21, F32, and points 1,2,3 in 3 images
		//without point normalization to be updated
		Trifocal(  Matr& Ma , Vectr& Ve , const Fundamental<T>& F21 ,const Fundamental<T>& F32, const std::vector<FVector<T,6> >& points ,T& beta ,T& gamma1,T& gamma2,T& gamma3)
		{

			Matrix3x3 f1 = F21.matrix();//.transpose();
			Matrix3x3 f2 = F32.matrix();//.transpose();

			Vector3 F21_1,F21_2,F21_3,F32_1 , F32_2,F32_3;
			F32_1[0] = f2(0,0);
			F32_1[1] = f2(1,0);
			F32_1[2] = f2(2,0);

			F32_2[0] = f2(0,1);
			F32_2[1] = f2(1,1);
			F32_2[2] = f2(2,1);

			F32_3[0] = f2(0,2);
			F32_3[1] = f2(1,2);
			F32_3[2] = f2(2,2);

			F21_1[0] = f1(0,0);
			F21_1[1] = f1(0,1);
			F21_1[2] = f1(0,2);

			F21_2[0] = f1(1,0);
			F21_2[1] = f1(1,1);
			F21_2[2] = f1(1,2);

			F21_3[0] = f1(2,0);
			F21_3[1] = f1(2,1);
			F21_3[2] = f1(2,2);

			//Calc epipoles normalisated
			Vector3 e1 = F21.rightEpipole();
			Vector3 e2 = F32.leftEpipole();

			//Calc basis change matrix
			Vector3 a2_1 = F21_1^e1;


			Vector3 a2_2 = F21_2^e1;
			Vector3 a2_3 = F21_3^e1;

			Vector3 b2_1 = F32_1^e2;
			Vector3 b2_2 = F32_2^e2;
			Vector3 b2_3 = F32_3^e2;

			Matrix3x3 A1,A2,A3,B1,B2,B3;

			A1(0,0) = F21_1[0];
			A1(0,1) = F21_1[1];
			A1(0,2) = F21_1[2];

			A1(1,0) = a2_1[0];
			A1(1,1) = a2_1[1];
			A1(1,2) = a2_1[2];

			A1(2,0) = e1[0];
			A1(2,1) = e1[1];
			A1(2,2) = e1[2];

			A2(0,0) = F21_2[0];
			A2(0,1) = F21_2[1];
			A2(0,2) = F21_2[2];

			A2(1,0) = a2_2[0];
			A2(1,1) = a2_2[1];
			A2(1,2) = a2_2[2];

			A2(2,0) = e1[0];
			A2(2,1) = e1[1];
			A2(2,2) = e1[2];

			A3(0,0) = F21_3[0];
			A3(0,1) = F21_3[1];
			A3(0,2) = F21_3[2];

			A3(1,0) = a2_3[0];
			A3(1,1) = a2_3[1];
			A3(1,2) = a2_3[2];

			A3(2,0) = e1[0];
			A3(2,1) = e1[1];
			A3(2,2) = e1[2];

			B1(0,0) = F32_1[0];
			B1(0,1) = F32_1[1];
			B1(0,2) = F32_1[2];

			B1(1,0) = b2_1[0];
			B1(1,1) = b2_1[1];
			B1(1,2) = b2_1[2];

			B1(2,0) = e2[0];
			B1(2,1) = e2[1];
			B1(2,2) = e2[2];

			B2(0,0) = F32_2[0];
			B2(0,1) = F32_2[1];
			B2(0,2) = F32_2[2];

			B2(1,0) = b2_2[0];
			B2(1,1) = b2_2[1];
			B2(1,2) = b2_2[2];

			B2(2,0) = e2[0];
			B2(2,1) = e2[1];
			B2(2,2) = e2[2];

			B3(0,0) = F32_3[0];
			B3(0,1) = F32_3[1];
			B3(0,2) = F32_3[2];

			B3(1,0) = b2_3[0];
			B3(1,1) = b2_3[1];
			B3(1,2) = b2_3[2];

			B3(2,0) = e2[0];
			B3(2,1) = e2[1];
			B3(2,2) = e2[2];
			
			//Calc normalisated list lines l2 , l3 with previous matrix
			Matr A(points.size() * 4,4);
			Vector<T> C(points.size()* 4);

			for (int i = 0; i<points.size(); i++)
			{	
				Vector3 u2,u3,v2,v3,p1;
				p1[0] = points[i].at(2);
				p1[1] = points[i].at(3);
				p1[2] = 1.0;
				u2[0] = 1.0;
				u3[0] = 1.0;
				u2[1] =0.0;
				u3[1]=0.0;
				u2[2]=- points[i].at(0);
				u3[2]=-points[i].at(4);
				v2[2]=-points[i].at(1);
				v3[2]=-points[i].at(5);
				v2[0]=0.0;
				v3[0]=0.0;
				v2[1]=1.0;
				v3[1]=1.0;
				//apply line according to change basis
				Vector3 u2_1,u3_1,v2_1,v3_1,u2_2,u3_2,v2_2,v3_2,u2_3,u3_3,v2_3,v3_3;
				u2_1 = A1 * u2;
				u2_2 = A2 * u2;
				u2_3 = A3 * u2;
				v2_1 = A1 * v2;
				v2_2 = A2 * v2;
				v2_3 = A3 * v2;

				u3_1 = B1 * u3;
				u3_2 = B2 * u3;
				u3_3 = B3 * u3;
				v3_1 = B1 * v3;
				v3_2 = B2 * v3;
				v3_3 = B3 * v3;

				//Compute matrix
				A(i*4,0) = p1[0] * u2_1[2] * u3_1[1] +  p1[1] * u2_2[2] * u3_2[1] +  p1[2] * u2_3[2] * u3_3[1];
				A(i*4,1) = p1[0] * u2_1[2] * u3_1[2];
				A(i*4,2) = p1[1] * u2_2[2] * u3_2[2];
				A(i*4,3) = p1[2] * u2_3[2] * u3_3[2];
				C(i*4) = -p1[0] * u2_1[1] * u3_1[2] -  p1[1] * u2_2[1] * u3_2[2] -  p1[2] * u2_3[1] * u3_3[2];

				A(i*4+1,0) = p1[0] * u2_1[2] * v3_1[1] +  p1[1] * u2_2[2] * v3_2[1] +  p1[2] * u2_3[2] * v3_3[1];
				A(i*4+1,1) = p1[0] * u2_1[2] * v3_1[2];
				A(i*4+1,2) = p1[1] * u2_2[2] * v3_2[2];
				A(i*4+1,3) = p1[2] * u2_3[2] * v3_3[2];
				C(i*4+1) = -p1[0] * u2_1[1] * v3_1[2] -  p1[1] * u2_2[1] * v3_2[2] -  p1[2] * u2_3[1] * v3_3[2];

				A(i*4+2,0) = p1[0] * v2_1[2] * u3_1[1] +  p1[1] * v2_2[2] * u3_2[1] +  p1[2] * v2_3[2] * u3_3[1];
				A(i*4+2,1) = p1[0] * v2_1[2] * u3_1[2];
				A(i*4+2,2) = p1[1] * v2_2[2] * u3_2[2];
				A(i*4+2,3) = p1[2] * v2_3[2] * u3_3[2];
				C(i*4+2) = -p1[0] * v2_1[1] * u3_1[2] - p1[1] * v2_2[1] * u3_2[2] -  p1[2] * v2_3[1] * u3_3[2];

				A(i*4+3,0) = p1[0] * v2_1[2] * v3_1[1] +  p1[1] * v2_2[2] * v3_2[1] +  p1[2] * v2_3[2] * v3_3[1];
				A(i*4+3,1) = p1[0] * v2_1[2] * v3_1[2];
				A(i*4+3,2) = p1[1] * v2_2[2] * v3_2[2];
				A(i*4+3,3) = p1[2] * v2_3[2] * v3_3[2];
				C(i*4+3) = -p1[0] * v2_1[1] * v3_1[2] -  p1[1] * v2_2[1] * v3_2[2] -  p1[2] * v2_3[1] * v3_3[2];
			}
			Ma = A;
			Ve = C;
			C = C / A.frobenius_norm();
			A = A / A.frobenius_norm();
			//Calc tenseur T_213 from fundamental matrices and trifocal points
			//solve b for A * b = C;
			Vector<T> b = A.linSolve(C);//give beta, gamma1,gamma2,gamma3.

			//update beta: useful to recover all cameras in one unique projective space
			beta = b[0];
			gamma1 = b[1];
			gamma2 = b[2];
			gamma3 = b[3];

			//Compute cameras
			FArray<T,3> bc;
			FMatrix<T,3,3> Ac;
			//tenseur 213
			//camera2 identity
			Ac.identity();
			bc.fill(T(0));
			_camera2 = Camera<T>( Ac,  bc);
			
			Vector3 t1 = a2_1 + b[1] * e1;
			Vector3 t2 = a2_2 + b[2] * e1;
			Vector3 t3 = a2_3 + b(3) * e1;
			Ac(0,0) = t1[0];
			Ac(1,0) = t1[1];
			Ac(2,0) = t1[2];

			Ac(0,1) = t2[0];
			Ac(1,1) = t2[1];
			Ac(2,1) = t2[2];

			Ac(0,2) = t3[0];
			Ac(1,2) = t3[1];
			Ac(2,2) = t3[2];

			_camera1 = Camera<T>(Ac , e1);

			//camera3
			Ac(0,0) = -b2_1[0]*b[0];
			Ac(1,0) = -b2_1[1]*b[0];
			Ac(2,0) = -b2_1[2]*b[0];

			Ac(0,1) = -b2_2[0]*b[0];
			Ac(1,1) = -b2_2[1]*b[0];
			Ac(2,1) = -b2_2[2]*b[0];

			Ac(0,2) = -b2_3[0]*b[0];
			Ac(1,2) = -b2_3[1]*b[0];
			Ac(2,2) = -b2_3[2]*b[0];

			_camera3 = Camera<T>(Ac , e2);

		}


		// Constructeur par recopie
		Trifocal& operator = (const Trifocal& trif) 
		{			
			_camera1 = trif._camera1;
			_camera2 = trif._camera2;
			_camera3 = trif._camera3;

			return *this;
		}


		// Access to matrix, camera, epipole
		const Camera<T> camera1() const { return _camera1; }
		const Camera<T> camera2() const { return _camera2; }
		const Camera<T> camera3() const { return _camera3; }


		//set parameters
		void setCamera1(Matrix3x3& A, FArray<T,3>& b)
		{
			_camera1 = Cam( A,  b);
		}
		void setCamera2(Matrix3x3& A, FArray<T,3>& b) 
		{
			_camera2 = Cam( A,  b);
		}
		void setCamera3(Matrix3x3& A, FArray<T,3>& b) 
		{
			_camera3 = Cam( A,  b);
		}




	protected:
		void computeCameras(Matr Teb1, Matr Teb2, Matr Teb3)
		{
			Vector<T> s,e1(3),e2(3),Ui,Vi;
			Matrix<T> U0(3,3),V0(3,3),U,Vt,V;
			FArray<T,3> b;
			FMatrix<T,3,3> Ac,B;
			FVector<T,3> K1,K2,K3;
			FMatrix<T,3,1> e2v;


			svd(Teb1,U,s,Vt);
			Vi=Vt.getRow(2);
			Ui=U.getCol(2);
			V0.setCol(0,Vi);
			U0.setCol(0,Ui);

			svd(Teb2,U,s,Vt);
			Vi=Vt.getRow(2);
			Ui=U.getCol(2);
			V0.setCol(1,Vi);
			U0.setCol(1,Ui);

			svd(Teb3,U,s,Vt);
			Vi=Vt.getRow(2);
			Ui=U.getCol(2);
			V0.setCol(2,Vi);
			U0.setCol(2,Ui);

			svd(U0,U,s,Vt);		
			e1 = U.getCol(2);
			//e1 = e1 / e1[2];
			svd(V0,U,s,Vt);
			e2 = U.getCol(2);
			//e2 = e2 / e2[2];

			Ac=FMatrix<T,3,3>::Identity();
			b.fill(T(0.));
			_camera1 = Camera<T>( Ac,  b);


			K1[0]=(Teb1*e2)[0];
			K1[1]=(Teb1*e2)[1];
			K1[2]=(Teb1*e2)[2];

			K2[0]=(Teb2*e2)[0];
			K2[1]=(Teb2*e2)[1];
			K2[2]=(Teb2*e2)[2];

			K3[0]=(Teb3*e2)[0];
			K3[1]=(Teb3*e2)[1];
			K3[2]=(Teb3*e2)[2];
			Ac.setCol(0,K1);
			Ac.setCol(1,K2);
			Ac.setCol(2,K3);

			b[0]=e1[0];
			b[1]=e1[1];
			b[2]=e1[2];

			_camera2 = Camera<T>( Ac,  b);
			Teb1 = transpose(Teb1);
			Teb2 = transpose(Teb2);
			Teb3 = transpose(Teb3);
			K1[0]=(Teb1*e1)[0];
			K1[1]=(Teb1*e1)[1];
			K1[2]=(Teb1*e1)[2];

			K2[0]=(Teb2*e1)[0];
			K2[1]=(Teb2*e1)[1];
			K2[2]=(Teb2*e1)[2];

			K3[0]=(Teb3*e1)[0];
			K3[1]=(Teb3*e1)[1];
			K3[2]=(Teb3*e1)[2];
			Ac.setCol(0,K1);
			Ac.setCol(1,K2);
			Ac.setCol(2,K3);


			e2v(0,0)=e2[0];
			e2v(1,0)=e2[1];
			e2v(2,0)=e2[2];


			b[0]=e2[0];
			b[1]=e2[1];
			b[2]=e2[2];
			B=FMatrix<T,3,3>::Identity();
			Ac = (e2v * transpose(e2v) - B) * Ac;
			_camera3 = Camera<T>( Ac,  b);
		}


		void computeMatrix(Matr &T1 , Matr &T2, Matr &T3 , Matr &X1_pts,Matr &X2_pts,Matr &X3_pts)
		{
			int n =0,r=15;
			Vector<T> s,Tensv(27),e1(3),e2(3),Ui,Vi;
			Matrix<T> U0(3,3),V0(3,3),U,Vt,V,Te1(3,3),Te2(3,3),Te3(3,3),Teb1(3,3),Teb2(3,3),Teb3(3,3),E(27,18),Ue(27,r),A(nbPoints*4,27), AU(nbPoints*4,r);
			FArray<T,3> b;
			FMatrix<T,3,3> Ac,B;
			FVector<T,3> K1,K2,K3;
			FMatrix<T,3,1> e2v;

			A.fill(0.);

			for(int i=0;i<nbPoints;i++)
			{
				n=4*i;
				A(n,0)=X1_pts(0,i);
				A(n,2)=- X1_pts(0,i) * X3_pts(0,i) ;
				A(n,6)=- X1_pts(0,i) * X2_pts(0,i) ;
				A(n,8)= X1_pts(0,i) * X2_pts(0,i) * X3_pts(0,i);
				A(n,9)= X1_pts(1,i);
				A(n,11)=- X1_pts(1,i) * X3_pts(0,i) ;
				A(n,15)=- X1_pts(1,i) * X2_pts(0,i) ;
				A(n,17)= X1_pts(1,i) * X2_pts(0,i) * X3_pts(0,i);
				A(n,18)= 1;
				A(n,20)= - X3_pts(0,i);
				A(n,24)= - X2_pts(0,i);
				A(n,26)=  X2_pts(0,i) * X3_pts(0,i);

				A(n+1,3)=X1_pts(0,i);
				A(n+1,5)=- X1_pts(0,i) * X3_pts(0,i) ;
				A(n+1,6)=- X1_pts(0,i) * X2_pts(1,i) ;
				A(n+1,8)= X1_pts(0,i) * X2_pts(1,i) * X3_pts(0,i);
				A(n+1,12)= X1_pts(1,i);
				A(n+1,14)=- X1_pts(1,i) * X3_pts(0,i) ;
				A(n+1,15)=- X1_pts(1,i) * X2_pts(1,i) ;
				A(n+1,17)= X1_pts(1,i) * X2_pts(1,i) * X3_pts(0,i);
				A(n+1,21)= 1;
				A(n+1,23)= - X3_pts(0,i);
				A(n+1,24)= - X2_pts(1,i);
				A(n+1,26)=  X2_pts(1,i) * X3_pts(0,i);

				A(n+2,1)=X1_pts(0,i);
				A(n+2,2)=- X1_pts(0,i) * X3_pts(1,i) ;
				A(n+2,7)=- X1_pts(0,i) * X2_pts(0,i) ;
				A(n+2,8)= X1_pts(0,i) * X2_pts(0,i) * X3_pts(1,i);
				A(n+2,10)= X1_pts(1,i);
				A(n+2,11)=- X1_pts(1,i) * X3_pts(1,i) ;
				A(n+2,16)=- X1_pts(1,i) * X2_pts(0,i) ;
				A(n+2,17)= X1_pts(1,i) * X2_pts(0,i) * X3_pts(1,i);
				A(n+2,19)= 1;
				A(n+2,20)= - X3_pts(1,i);
				A(n+2,25)= - X2_pts(0,i);
				A(n+2,26)=  X2_pts(0,i) * X3_pts(1,i);

				A(n+3,4)=X1_pts(0,i);
				A(n+3,5)=- X1_pts(0,i) * X3_pts(1,i) ;
				A(n+3,7)=- X1_pts(0,i) * X2_pts(1,i) ;
				A(n+3,8)= X1_pts(0,i) * X2_pts(1,i) * X3_pts(1,i);
				A(n+3,13)= X1_pts(1,i);
				A(n+3,14)=- X1_pts(1,i) * X3_pts(1,i) ;
				A(n+3,16)=- X1_pts(1,i) * X2_pts(1,i) ;
				A(n+2,17)= X1_pts(1,i) * X2_pts(1,i) * X3_pts(1,i);
				A(n+3,22)= 1;
				A(n+3,23)= - X3_pts(1,i);
				A(n+3,25)= - X2_pts(1,i);
				A(n+3,26)=  X2_pts(1,i) * X3_pts(1,i);

			}
			svd(A,U,s,Vt,true);
			Tensv = Vt.getRow(26);
			/*Tenseur Trifocal-----------------------------
			----------------------------------------------*/
			Te1(0,0)=Tensv[0];
			Te1(0,1)=Tensv[1];
			Te1(0,2)=Tensv[2];
			Te1(1,0)=Tensv[3];
			Te1(1,1)=Tensv[4];
			Te1(1,2)=Tensv[5];
			Te1(2,0)=Tensv[6];
			Te1(2,1)=Tensv[7];
			Te1(2,2)=Tensv[8];

			Te2(0,0)=Tensv[9];
			Te2(0,1)=Tensv[10];
			Te2(0,2)=Tensv[11];
			Te2(1,0)=Tensv[12];
			Te2(1,1)=Tensv[13];
			Te2(1,2)=Tensv[14];
			Te2(2,0)=Tensv[15];
			Te2(2,1)=Tensv[16];
			Te2(2,2)=Tensv[17];

			Te3(0,0)=Tensv[18];
			Te3(0,1)=Tensv[19];
			Te3(0,2)=Tensv[20];
			Te3(1,0)=Tensv[21];
			Te3(1,1)=Tensv[22];
			Te3(1,2)=Tensv[23];
			Te3(2,0)=Tensv[24];
			Te3(2,1)=Tensv[25];
			Te3(2,2)=Tensv[26];

			/*Find Epipoles-----------------------------
			----------------------------------------------*/
			svd(Te1,U,s,Vt);
			Vi=Vt.getRow(2);
			Ui=U.getCol(2);
			V0.setCol(0,Vi);
			U0.setCol(0,Ui);

			svd(Te2,U,s,Vt);
			Vi=Vt.getRow(2);
			Ui=U.getCol(2);
			V0.setCol(1,Vi);
			U0.setCol(1,Ui);

			svd(Te3,U,s,Vt);
			Vi=Vt.getRow(2);
			Ui=U.getCol(2);
			V0.setCol(2,Vi);
			U0.setCol(2,Ui);
			svd(U0,U,s,Vt);		
			e1 = U.getCol(2);
			//e1 = e1 / e1[2];
			svd(V0,U,s,Vt);
			e2 = U.getCol(2);
			//e2 = e2 / e2[2];

			/*Run Minimization from Epipoles-----------------------------
			----------------------------------------------*/
			E.fill(0.);

			for(int i=0;i<3;i++)
			{
				E(0+9*i,0+6*i)=e2[0];
				E(1+9*i,0+6*i)=e2[1];
				E(2+9*i,0+6*i)=e2[2];	

				E(3+9*i,1+6*i)=e2[0];
				E(4+9*i,1+6*i)=e2[1];
				E(5+9*i,1+6*i)=e2[2];	

				E(6+9*i,2+6*i)=e2[0];
				E(7+9*i,2+6*i)=e2[1];
				E(8+9*i,2+6*i)=e2[2];	

				E(0+9*i,3+6*i)=-e1[0];
				E(1+9*i,4+6*i)=-e1[0];
				E(2+9*i,5+6*i)=-e1[0];

				E(3+9*i,3+6*i)=-e1[1];
				E(4+9*i,4+6*i)=-e1[1];
				E(5+9*i,5+6*i)=-e1[1];

				E(6+9*i,3+6*i)=-e1[2];
				E(7+9*i,4+6*i)=-e1[2];
				E(8+9*i,5+6*i)=-e1[2];
			}
			svd(E,U,s,Vt,true);	
			//cout << "U " << U << endl;
			for (int i=0;i<r;i++)
			{
				Ue.setCol(i,U.getCol(i));
			}
			AU = A * Ue;
			svd(AU,U,s,Vt,true);
			Tensv = Ue *  Vt.getRow(r-1);

			/*Compute again Tenseur-----------------------------
			----------------------------------------------*/
			Te1(0,0)=Tensv[0];
			Te1(0,1)=Tensv[1];
			Te1(0,2)=Tensv[2];
			Te1(1,0)=Tensv[3];
			Te1(1,1)=Tensv[4];
			Te1(1,2)=Tensv[5];
			Te1(2,0)=Tensv[6];
			Te1(2,1)=Tensv[7];
			Te1(2,2)=Tensv[8];

			Te2(0,0)=Tensv[9];
			Te2(0,1)=Tensv[10];
			Te2(0,2)=Tensv[11];
			Te2(1,0)=Tensv[12];
			Te2(1,1)=Tensv[13];
			Te2(1,2)=Tensv[14];
			Te2(2,0)=Tensv[15];
			Te2(2,1)=Tensv[16];
			Te2(2,2)=Tensv[17];

			Te3(0,0)=Tensv[18];
			Te3(0,1)=Tensv[19];
			Te3(0,2)=Tensv[20];
			Te3(1,0)=Tensv[21];
			Te3(1,1)=Tensv[22];
			Te3(1,2)=Tensv[23];
			Te3(2,0)=Tensv[24];
			Te3(2,1)=Tensv[25];
			Te3(2,2)=Tensv[26];
			/*Unormalize-----------------------------
			----------------------------------------------*/
			T2 = inverse(T2);
			T3 = transpose(inverse(T3));

			Te1 = T2 * Te1 * T3;
			Te2 = T2 * Te2 * T3;
			Te3 = T2 * Te3 * T3;

			Teb1 =  Te1 * T1(0,0) + Te2 * T1(1,0) + Te3 * T1(2,0);
			Teb2 =  Te1 * T1(0,1) + Te2 * T1(1,1) + Te3 * T1(2,1);
			Teb3 =  Te1 * T1(0,2) + Te2 * T1(1,2) + Te3 * T1(2,2);
			/*Find Epipole for Unormalized and minimized tensor-----
			-------------------------------------------------------*/
			computeCameras(Teb1, Teb2, Teb3);


		}


		Matr normalize2D(Matr &X1_pts) 
		{
			Matr X1_ptsn(3,nbPoints)  ;
			Matr T1=Matr::Zero(3,3) ;
			Vectr c_mean=Vectr::Zero(3);
			int i;
			T scale=0;


			for (i=0; i<nbPoints; i++)
			{
				c_mean +=  X1_pts.getCol(i);
			}
			c_mean = c_mean / T(nbPoints);
			for (i=0; i<nbPoints; i++)
			{
				X1_ptsn(0,i) = X1_pts(0,i) - c_mean[0];
				X1_ptsn(1,i) = X1_pts(1,i) - c_mean[1];
			}

			for (i=0; i<nbPoints; i++)
			{
				scale +=  std::sqrt(pow(X1_ptsn(0,i),2) + pow(X1_ptsn(1,i),2));
			}
			scale = scale  / T(nbPoints);
			scale = std::sqrt((T)2) / scale;

			T1(0,0)=scale;
			T1(1,1)=scale;
			T1(2,2)=1;
			T1(0,2)= - scale * c_mean[0];
			T1(1,2)= - scale * c_mean[1];

			return T1;
		}


	};
	/// @}
}

