// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{


	// Functor: estimate line equation: Tensor matrix from a sample of 7 points in 3 images
	template <typename V ,class Point_Format>
	class TrifEstimator {
	public:
		template <class OutputTensorIterator>
		void operator() (const FArray<Point_Format,7>& points,
			OutputTensorIterator lines) const {
				Trifocal<V,7,Point_Format> line(points);
				*lines++ = line;
		}
	};

	// Functor: estimate the square residual of a point relatively to a line equation
	template <typename V ,class Point_Format>
	class TrifResidual {
	public:
		template <class T>
		V operator() (const Trifocal<V,7,Point_Format>& line,
			const T& point) const {
				Camera<V> P1 = line.camera1();
				Camera<V> P2 = line.camera2();
				Camera<V> P3 = line.camera3();

				FVector<V,2> p1,p2,p3;

				p1[0]=point[0];
				p1[1]=point[1];

				p2[0]=point[2];
				p2[1]=point[3];

				p3[0]=point[4];
				p3[1]=point[5];
				FVector<V,3> Pt;

				Triangulation<V> triangulation;
				triangulation.add(P1,p1);
				triangulation.add(P2,p2);
				triangulation.add(P3,p3);
				triangulation.compute();


				return triangulation.error();
		}
	};
	/// @}
}
