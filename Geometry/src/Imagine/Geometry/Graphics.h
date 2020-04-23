// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{

	template <typename T>
	void getCamera(Camera<T>& camera, int width, int height) {

		FVector<T,3> pos,dir,up;
		T zoom;
		getCamera(pos,dir,up,zoom);
		dir.normalize();
		up.normalize();

		typename Camera<T>::Matrix3x3 matrix;
		matrix.setRow(2,dir);
		matrix.setRow(1,-up);
		matrix.setRow(0,-up^dir);

		T pixelSize = T((height-1) / 2) / tan(T(M_PI * 30 / 360) / zoom);

		camera = Camera<T>(
			Affine<T,3>( matrix, - matrix * pos ),
			T((width-1) / 2),
			T((height-1) / 2),
			pixelSize,
			pixelSize);
	}
	/// @}
}
