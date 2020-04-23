// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

#include <fstream>
#include <iostream>
#include <vector>

#include <Imagine/Features.h>

using namespace std;
using namespace Imagine;

template <typename T>
void detection(const T& detector) {
	cout << "Detection ... ";
	Image<byte> I;
	load(I,srcPath("toys.0000.png"));
	Array<typename T::Feature> feats=detector.run(I);	// runs detection
	cout << feats.size() << " points" << endl;

}	

// Good/bad criterion based on fundamental matrix.
class FCrit {
	FMatrix<float,3,3> F;
	float thres;
public:
	FCrit(const FMatrix<float,3,3>& f, float th): F(f), thres(th) {}
	inline bool operator()(const Array<SIFTDetector::Feature>& F1,size_t i1,const Array<SIFTDetector::Feature>& F2,size_t i2) const {
		FloatPoint2 a=F1[i1].pos,b=F2[i2].pos;
		FVector<float,3> m1(a.x(), a.y(), 1.0f), m2(b.x(), b.y(), 1.0f);
		FVector<float,3> line = F*m1;
		float d1 = abs((line*m2)/sqrt(line[0]*line[0]+line[1]*line[1]));
		line = tmult(F,m2);
		float d2 = abs((line*m1)/sqrt(line[0]*line[0]+line[1]*line[1]));
		return max(d1,d2)<thres;
	}
};
// ...

void matching() {
	cout << "Matching" << endl;
	Image<byte> I1,I2;
	load(I1,srcPath("toys.0000.png"));
	load(I2,srcPath("toys.0200.png"));
	
	SIFTDetector d;
	Array<SIFTDetector::Feature> feats1=d.run(I1);	
	Array<SIFTDetector::Feature> feats2=d.run(I2);	
	ifstream str(srcPath("toys_F_0200_0000.txt"));
	FMatrix<double,3,3> F;
	str >> F;
	F=transpose(F); // From 0 to 200
	str.close();
	
	size_t bi2;
	double bd,bd2;
	find2Best(feats1,feats2,0,bi2,bd,bd2);						// two best matches
	cout << "0 -> " << bi2 << " ratio= " << bd/bd2 << endl;
	find2Best(feats1,feats2,0,bi2,bd,bd2,&F,10.);		
	cout << "0 -> " << bi2 << " ratio= " << bd/bd2 << endl;		// ...

	MatchList L1,L2;									// Lowe's match
	L1=loweMatch(feats1,feats2,.6);
	L2=loweMatch(feats1,feats2,.75,true,200,&F,5.);
	cout << L1.size() << " " << L2.size() << endl;		// ...

	int w=I1.width(),h=I1.height();
	openWindow(w,2*h);
	display(I1);display(I2,0,h);
	writeFeaturePoints(feats1,"keys1.txt");			// write FP
	writeFeaturePoints(feats2,"keys2.txt",true);	// ...
	readFeaturePoints(feats1,"keys1.txt");			// read FP
	readFeaturePoints(feats2,"keys2.txt",true);		// ...
	
	drawFeatures(feats1);													// draw features
	drawFeatures(feats2,IntPoint2(0,h),BLUE);									// ...
	drawFeature(feats1[10],IntPoint2(0,0),GREEN,true,true);					// draw one feature
	click();
	display(I1);display(I2,0,h);
	drawMatches(feats1,feats2,L2,IntPoint2(0,0),IntPoint2(0,h),1.f,true);		// draw matches
	click();
	display(I1);display(I2,0,h);
	drawGoodMatches(feats1,feats2,L1,FCrit(F,1.f),IntPoint2(0,0),IntPoint2(0,h),1.f,true);	// draw good matches
	click();
}	

int main() {
	SIFTDetector sd;			// sift, constructor
	sd.setNumOctaves(5);		// sift, number of octaves
	sd.setFirstOctave(0);		// sift, first octave
	sd.setNumScales(4);			// sift, number of scales per octave
	sd.setEdgeThresh(20.f);		// sift, max ratio of Hessian eigenvalues
	sd.setPeakThresh(0.05f);	// sift, min contrast
	detection(sd);
	
	LoweSIFTDetector lsd;		// Lowe sift, constructor
	lsd.setNumOctaves(5);		// Lowe sift, number of octaves
	lsd.setFirstOctave(0);		// Lowe sift, first octave
	lsd.setNumScales(4);		// Lowe sift, number of scales per octave
	lsd.setEdgeThresh(20.f);	// Lowe sift, max ratio of Hessian eigenvalues
	lsd.setPeakThresh(0.05f);	// Lowe sift, min contrast
	detection(lsd);
	
	MSERDetector md;			// MSER, constructor
	detection(md);

	matching();

	
	
	//  endGraphics();
    return 0;
}
