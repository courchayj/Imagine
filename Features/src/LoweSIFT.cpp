#include <iostream>
#include <list>

#include "Imagine/Features.h"

#include "lowe/LoweUtil.h"

using namespace std;
using namespace Imagine;

namespace Imagine {

	Array<LoweSIFT> LoweSIFTDetector::run(const Image<byte>& I) const {
		int w = I.width(), h = I.height(); 
		// initialize the image
		LoweImage lowe_img = CreateLoweImage(h, w, IMAGE_POOL); 
		for (int r = 0; r < h; r++) 
			for (int c = 0; c < w; c++)
				lowe_img->pixels[r][c] = I(c,r)/255.f;
		// compute the key points
		LoweKeypoint lowe_keys =
            GetLoweKeypoints(lowe_img,
                             firstOctave,numOctaves,numScales,
                             edgeThresh,peakThresh);
		list<LoweSIFT> L;
		while(lowe_keys != NULL)
		{
			LoweSIFT m;
			m.pos=FloatPoint2(lowe_keys->col,lowe_keys->row);
			m.scale=lowe_keys->scale;
			m.angle=lowe_keys->ori;
			m.desc.copy(lowe_keys->ivec);
			L.push_back(m);
			lowe_keys = lowe_keys->next;
		}

		// free memory
		FreeStoragePool(IMAGE_POOL);
		FreeStoragePool(KEY_POOL);

		return(Array<LoweSIFT>(L));
	}
}
