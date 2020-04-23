// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine{ 
	/// \addtogroup Geometry
	/// @{


/// Normalize set of points in dimension \a dim so that their center is 0 and
/// variance is 2. Return the applied similarity matrix.
template <typename T, int dim, typename Iterator>
FMatrix<T,dim+1,dim+1> normalize(Iterator begin, Iterator end)
{
    Iterator it=begin;
    // Compute mean
    FVector<T,dim> mean(T(0));
    int n=0;
    for(; it != end; ++it, ++n)
        mean += *it;
    mean /= T(n);
    // Center points and compute variance
    T var=T(0);
    for(it=begin; it != end; ++it) {
        *it -= mean;
        var += norm2(*it);
    }
    // Set variance=2
	var = std::sqrt( T(2*n)/var );
    for(it=begin; it != end; ++it)
        *it *= var;
    // Fill normalization matrix
    FMatrix<T,dim+1,dim+1> N(T(0));
    for(int i=0; i < dim; i++) {
        N(i,i) = var;
        N(i,dim) = -mean[i]*var;
    }
    N(dim,dim) = T(1);
    return N;
}

	/// @}
}

