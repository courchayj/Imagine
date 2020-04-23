// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

#include <vector>

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Camera triangulation using the iterated linear method

	template <typename T>
	class Triangulation
	{
	public:
		typedef Camera<T> CameraT;
		typedef FVector<T,2> Vector2;
		typedef FVector<T,3> Vector3;
		typedef FMatrix<T,3,3> Matrix3x3;


		void clear()
		{
			views.clear();
		}


		void add(const CameraT& camera, Vector2 p)
		{
			views.push_back( std::pair<const CameraT*, Vector2>(&camera,p) );
		}


		T error(Vector3 P) const
		{
			T squared_reproj_error = 0;
			for (std::size_t i=0;i<views.size();i++)
			{
				const CameraT& camera = *(views[i].first);
				Vector2 p = views[i].second;
				T x,y;
				camera.projection(P.x(),P.y(),P.z(),x,y);
				T dx = x-p.x(), dy = y-p.y();
				squared_reproj_error += dx*dx+dy*dy;
			}
			return squared_reproj_error;
		}


		Vector3 compute(int iter = 3) const
		{
			int nviews = int(views.size());
			assert(nviews>=2);

			// Iterative weighted linear least squares
			Matrix3x3 AtA;
			Vector3 Atb, P;
			std::vector<T> weights(nviews,T(1));
			for (int it=0;it<iter;it++)
			{
				AtA.fill(0);
				Atb.fill(0);
				for (int i=0;i<nviews;i++)
				{
					const CameraT& camera = *(views[i].first);
					const Vector2 p = views[i].second;
					const T w = weights[i];

					Vector3 v1,v2;
					for (int j=0;j<3;j++)
					{
						v1[j] = w * ( camera.matrix()(0,j) - p.x() * camera.matrix()(2,j) );
						v2[j] = w * ( camera.matrix()(1,j) - p.y() * camera.matrix()(2,j) );
						Atb[j] += w * ( v1[j] * ( p.x() * camera.vector()[2] - camera.vector()[0] ) + v2[j] * ( p.y() * camera.vector()[2] - camera.vector()[1] ) );
					}

					for (int k=0;k<3;k++) 
					{
						for (int j=0;j<=k;j++)
						{
							T v = v1[j] * v1[k] + v2[j] * v2[k];
							AtA(j,k) += v;
							if (j<k) AtA(k,j) += v;
						}
					}
				}

				P = inverse(AtA) * Atb;

				// Compute reprojection error, min and max depth, and update weights
				zmin = std::numeric_limits<T>::max();
				zmax = - std::numeric_limits<T>::max();
				err = 0;
				for (int i=0;i<nviews;i++)
				{
					const CameraT& camera = *(views[i].first);
					Vector2 p = views[i].second;
					T x,y,z;
					camera.projection(P.x(),P.y(),P.z(),x,y,z);
					if (z<zmin) zmin=z;
					if (z>zmax) zmax=z;
					T dx = x-p.x(), dy = y-p.y();
					err += dx*dx+dy*dy;
					weights[i] = T(1) / z;
				}
			}
			return P;
		}


		//Vector3 compute(T& squared_reproj_error, bool& in_front, T& rcond) const
		//{
		//	int nviews = int(views.size());
		//	if (nviews<2)
		//	{
		//		std::cerr << "Insufficient number of views for triangulation" << std::endl;
		//		return Vector3(0,0,0);
		//	}

		//	CL::LinAlg::Matrix<T> A(2*nviews,3);
		//	CL::LinAlg::Vector<T> b(2*nviews), weights(nviews);
		//	Vector3 P;

		//	// Iterative weighted least squares
		//	weights.fill(T(1));
		//	for (int it=0;it<iter;it++)
		//	{
		//		// Compute linear system coefficients
		//		for (int i=0;i<nviews;i++)
		//		{
		//			const Camera& camera = *(views[i].first);
		//			Vector2 p = views[i].second;
		//			const float w = weights(i);
		//			for (int j=0;j<3;j++)
		//			{
		//				A(2*i  ,j) = w * ( camera.matrix()(0,j) - p.x() * camera.matrix()(2,j) );
		//				A(2*i+1,j) = w * ( camera.matrix()(1,j) - p.y() * camera.matrix()(2,j) );
		//			}
		//			b(2*i  ) = w * ( p.x() * camera.vector()(2) - camera.vector()(0) );
		//			b(2*i+1) = w * ( p.y() * camera.vector()(2) - camera.vector()(1) );
		//		}

		//		// Solve linear system in least squares sense
		//		CL::LinAlg::Vector<T> P_ = A.linSolve(b);
		//		P = Vector3(P_(0),P_(1),P_(2));

		//		// Compute reprojection error and update weights
		//		squared_reproj_error = 0;
		//		in_front = true;
		//		for (int i=0;i<nviews;i++)
		//		{
		//			const Camera& camera = *(views[i].first);
		//			Vector2 p = views[i].second;
		//			T x,y,z;
		//			camera.projection(P.x(),P.y(),P.z(),x,y,z);
		//			if (z<=0) in_front = false;
		//			T dx = x-p.x(), dy = y-p.y();
		//			squared_reproj_error += dx*dx+dy*dy;
		//			weights(i) = T(1) / z;
		//		}
		//	}

		//	rcond = A.rcond();
		//	return P;
		//}


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Accessors

		// These values are defined after a successful call to compute
		T minDepth() const { return zmin; }
		T maxDepth() const { return zmax; }
		T error() const { return err; }

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Data members

	protected:
		mutable T zmin, zmax, err;
		std::vector< std::pair<const CameraT*, Vector2> > views;
	};
	/// @}
}


