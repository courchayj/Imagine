// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{

	// Imports various types of camera files
	template <typename T>
	bool readVistaCamera(Camera<T>& camera, std::istream &in)
	{
		// Lecture du fichier
		int nsels_x = 512, nsels_y = 480;
		int npixels_x = 512, npixels_y = 480;
		float sel_size_x = 0, sel_size_y = 0;
		float pixel_size_x = 0, pixel_size_y = 0;
		float image_center_x = 0, image_center_y = 0;
		float scale_x = 1;
		float rotation_x = 0, rotation_y = 0, rotation_z = 0;
		float translation_x = 0, translation_y = 0, translation_z = 0;
		float focal_length = 1;
		std::string s;
		do {
			getline(in,s);
#define READ_INT(key,val)   if (s.find(key)!=std::string::npos) val = atoi ( s.substr(s.find_first_of("-.0123456789")).c_str())
#define READ_FLOAT(key,val) if (s.find(key)!=std::string::npos) val = (float)atof ( s.substr(s.find_first_of("-.0123456789")).c_str())
			READ_INT("nsels_x",nsels_x);
			READ_INT("nsels_y",nsels_y);
			READ_INT("npixels_x",npixels_x);
			READ_INT("npixels_y",npixels_y);
			READ_FLOAT("sel_size_x",sel_size_x);
			READ_FLOAT("sel_size_y",sel_size_y);
			READ_FLOAT("pixel_size_x",pixel_size_x);
			READ_FLOAT("pixel_size_y",pixel_size_y);
			READ_FLOAT("image_center_x",image_center_x);
			READ_FLOAT("image_center_y",image_center_y);
			READ_FLOAT("scale_x",scale_x);
			READ_FLOAT("rotation_x",rotation_x);
			READ_FLOAT("rotation_y",rotation_y);
			READ_FLOAT("rotation_z",rotation_z);
			READ_FLOAT("translation_x",translation_x);
			READ_FLOAT("translation_y",translation_y);
			READ_FLOAT("translation_z",translation_z);
			READ_FLOAT("focal_length",focal_length);
#undef READ_INT
#undef READ_FLOAT
		} while (!in.eof());

		// Valeurs par defaut pour certains parametres
		if (pixel_size_x == 0) pixel_size_x = sel_size_x * nsels_x / npixels_x;
		if (pixel_size_y == 0) pixel_size_y = sel_size_y * nsels_y / npixels_y;
		if (image_center_x == 0) image_center_x = float(npixels_x) / 2;
		if (image_center_y == 0) image_center_y = float(npixels_y) / 2;

		// Construction de la matrice de projection
		camera = Camera<T>(
			Affine<T,3>::Translation(T(translation_x),T(translation_y),T(translation_z)) * Affine<T,3>::Rotation(T(rotation_x),T(rotation_y),T(rotation_z)),
			T(image_center_x),
			T(image_center_y),
			T(focal_length / pixel_size_x * scale_x),
			T(-focal_length / pixel_size_y)
			);

		return true;
	}


	template <typename T>
	bool readLfmCamera(Camera<T>& camera, std::istream &in)
	{
		T fov,nx,ny;
		in >> fov >> nx >> ny;
		typename Camera<T>::Matrix3x3 matrix;
		typename Camera<T>::Vector3 vector;
		in >> vector;
		in >> matrix;
		matrix.setRow(1, -matrix.getRow(1));
		matrix.setRow(2, -matrix.getRow(2));
		vector = - matrix * vector;
		T pixelSize = camera.deltav / tan(float(M_PI) * fov / 360);
		camera = Camera<T>(
			Affine<T,3>(matrix,vector),
			(nx-1) / 2,
			(ny-1) / 2,
			pixelSize,
			pixelSize
			);
		return true;
	}


	template <typename T>
	bool readKRTCamera(Camera<T>& camera, std::istream &in)
	{
		typename Camera<T>::Matrix3x3 K,R;
		typename Camera<T>::Vector3 t;
		in >> K;
		in >> R;
		in >> t;
		camera = Camera<T>(K * R, - K * R * t);
		return true;
	}


	template <typename T>
	bool readKULeuvenCamera(Camera<T>& camera, std::istream &in)
	{
		typename Camera<T>::Matrix3x3 K,R;
		typename Camera<T>::Vector3 t;
		in >> K;
		in >> t;
		in >> R;
		in >> t;
		camera = Camera<T>(K * ranspose(R), - K * transpose(R) * t);
		return true;
	}
	/// @}
}
