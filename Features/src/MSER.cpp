#include <iostream>
#include <list>

#include "Imagine/Features.h"

using namespace std;

extern "C" {
#include "vl/mser.h"
#include "vl/sift.h"
}

namespace Imagine {

	typedef FeaturePoint<byte> PNS; // Just pos and scale.
	bool scaleOrder(PNS m1,PNS m2) {
		return (m1.scale < m2.scale);
	}

	Array<MSER> MSERDetector::run(const Image<byte>& I) const {
		int w=I.width(),h=I.height();

		double   delta         = -1 ;
		double   max_area      = -1 ;
		double   min_area      = -1 ;
		double   max_variation = -1 ;
		double   min_diversity = -1 ;

		int dims[]={w,h};
		VlMserFilt *filt=vl_mser_new(2,dims);
		if (delta         >= 0) vl_mser_set_delta          (filt, (vl_mser_pix) delta) ;
		if (max_area      >= 0) vl_mser_set_max_area       (filt, max_area) ;
		if (min_area      >= 0) vl_mser_set_min_area       (filt, min_area) ;
		if (max_variation >= 0) vl_mser_set_max_variation  (filt, max_variation) ;
		if (min_diversity >= 0) vl_mser_set_min_diversity  (filt, min_diversity) ;

		vl_mser_process (filt,I.data()) ;

		filt->acc=0;	// So that ell_fit will alloc it
		vl_mser_ell_fit (filt) ;
		int n=vl_mser_get_ell_num(filt);
		int dof= vl_mser_get_ell_dof(filt);
		assert(dof==5);

		// 
		Array<PNS> pns(n);
		float const *frames=vl_mser_get_ell(filt);
		for (int i=0;i<n;++i,frames+=dof) {
			pns[i].pos=FloatPoint2(frames[0],frames[1]);
			// could be geometric mean of ell axis. Here sqrt(mean) seems more invariant?!
			pns[i].scale=max(1.f,3*std::sqrt(std::sqrt(std::sqrt(abs(frames[2]*frames[4]-frames[3]*frames[3]))))); 
		}
		vl_mser_delete (filt) ;
		// Now, use SIFT to compute orientation and descriptor
		Image<float,2> If(I);
		VlSiftFilt *sfilt=vl_sift_new (w,h,-1,3,0);
		bool first=true;

		sort(pns.data(),pns.data()+n,scaleOrder);
		int idx=0; // points are sorted.  
		list<MSER> L;
		while (true) {
			if (first) {
				first=false ;
				vl_sift_process_first_octave (sfilt, If.data());
			} else if (vl_sift_process_next_octave(sfilt))
				break; // Last octave
			int oct=vl_sift_get_octave_index (sfilt);
			for (;idx<n;idx++) {
				VlSiftKeypoint ik ;
				vl_sift_keypoint_init (sfilt,&ik,pns[idx].x(), pns[idx].y(),pns[idx].scale);
				if (ik.o > oct)
					break;	// Keep the the remaining ones for next octaves
				else if (ik.o < oct)
					continue;
				double angles[4]; // Use SIFT ori
				int  nangles = vl_sift_calc_keypoint_orientations(sfilt,angles,&ik) ;            
				MSER m;
				m.pos=pns[idx].pos;
				m.scale=pns[idx].scale;
				for (int q = 0 ; q < nangles ; ++q) {
					m.angle=float(angles[q]);			
					vl_sift_pix descr[128];
					vl_sift_calc_keypoint_descriptor(sfilt,descr,&ik,angles [q]) ;
					for (int k=0;k<128;k++)
						m.desc[k]=byte(512*descr[k]);
					L.push_back(m);
				}
			}
		}
		vl_sift_delete(sfilt);
		return Array<MSER>(L); 
	}


}

