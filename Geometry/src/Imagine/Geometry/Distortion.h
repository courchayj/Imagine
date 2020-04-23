// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{



	//	// ================================== formules en r^3 ===============================
	//	// Distortion
	//	// Taken from Reg Willson's implementation of Tsai's calibration technique.
	//	// The routine algebraically solves a cubic polynomial in Rd using the Cardan method.
	//#define SQRT3   1.732050807568877293527446341505872366943
	//#define CUB(X)((X)*(X)*(X))
	//#define SQR(X)((X)*(X))
	//#if WIN32
	//#define cbrt(x) (((x)>0)?pow((x),1/3.0):-pow(-(x),1/3.0))
	//#define hp(x,y) _hypot((x),(y))
	//#else
	//#define hp(x,y) hypot((x),(y))
	//#endif
	//
	//	void undis2disSensorCoords(double disto,double Xu,double Yu,double &Xd,double &Yd)  {
	//		double Ru, Rd, lambda, c, d, Q, R, D, S, T, sinT, cosT;
	//		// first check for boundary conditions
	//		if (((Xu == 0) && (Yu == 0)) || (disto == 0.0)) {
	//			Xd = Xu;
	//			Yd = Yu;
	//			return;
	//		}
	//		Ru = hp(Xu, Yu);        // sqrt(Xu*Xu+Yu*Yu)
	//		// find the "good" real root of the equation Ru = (1 + kappa1 * Rd**2)*Rd in Rd
	//		// the equation is x**3 + c*x + d = 0
	//		// http://lib-www.lanl.gov/numerical/bookcpdf/c5-6.pdf
	//		// Cubic Equations
	//		c = 1. / disto;
	//		d = -c * Ru;
	//		Q = c / 3.;
	//		R = -d / 2.;
	//		D = CUB (Q) + SQR (R);
	//		if (D >= 0) {               // in most cases, D >= 0
	//			D = sqrt (D);
	//			S = cbrt (R + D);
	//			T = cbrt (R - D);
	//			Rd = S + T;
	//		} else {                    // D < 0 => 3 real roots
	//			D = sqrt (-D);
	//			S = cbrt (hp(R, D));
	//			T = atan2 (D, R) / 3.;
	//			sinT = sin(T);
	//			cosT = cos(T);
	//			Rd = -S * cosT + SQRT3 * S * sinT;
	//			// other roots: 2.0*S*cos(T), -S*cos(T)-SQRT3*S*sin(T)
	//		}
	//		lambda = Rd / Ru;
	//		//* Horrible hack
	//		if(lambda < 0.0){
	//			std::cout << "oops disto" << std::endl;
	//			lambda = 1.0;
	//		}
	//		Xd = Xu * lambda;
	//		Yd = Yu * lambda;
	//	}


	// Main function: undistort image Id to Iu, according to (cam,dis)
	template <typename TI,typename TC>
	Image<TI> correctRadialDistortion(const Image<TI>& Id, const Camera<TC>& cam) {
		int w=Id.width(),h=Id.height();
		Image<TI> Iu(w,h);
		for (int y=0;y<h;y++)
			for (int x=0;x<w;x++) {
				TC xd,yd;
				cam.addDistortion(TC(x),TC(y),xd,yd);
				if (xd<0 || yd<0 || xd>w-1 || yd>h-1)
					Iu(x,y)=TI( typename PixelTraits<TI>::scalar_type(0) );
				else
					Iu(x,y)=TI( Id.interpolate(xd,yd) );
			}
			return Iu;
	}

	//// ========================= Formules a la PtGUI ========================

	//// Main function: undistort image Id to Iu, according to (xc,yc,a,b,c)
	//// Warning: ne marche que si TU repose sur un type flottant... A changer!!
	//template <typename TU,typename TD,typename T>
	//Images::Image<TU,2> correctRadialDistortion(const Images::Image<TD,2>& Id,T a,T b,T c,T xc=0,T yc=0) {
	//	int w=Id.width(),h=Id.height();
	//	T d=1-a-b-c;
	//	xc+=w/T(2);yc+=h/T(2);
	//	T n=min(w/T(2),h/T(2));
	//	Images::Image<TU,2> Iu(w,h);
	//	for (int y=0;y<h;y++)
	//		for (int x=0;x<w;x++) {
	//			T xd=x-xc,yd=y-yc;
	//			T r=sqrt(xd*xd+yd*yd)/n;
	//			T f=(a*r*r*r+b*r*r+c*r+d);
	//			xd=xc+f*xd;
	//			yd=yc+f*yd;
	//			if (xd<0 || yd<0 || xd>w-1 || yd>h-1)
	//				Iu(x,y)=TU(typename Images::Scalar_type<TU>::scalar_type(0));
	//			else
	//				Iu(x,y)=Images::interpol<TU>(Id,typename Images::Scalar_type<TU>::scalar_type(xd),typename Images::Scalar_type<TU>::scalar_type(yd));
	//		}
	//	return Iu;
	//}

	/// @}
}

