// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

// No arg: Example testing features
// 2 args: ima output (extracts from ima and writes SIFTs to output, Lowe format)
// 3 args: i1 output i2 (extracts from i1 and i2 and write matches to output)

#include <fstream>
#include <iostream>
#include <vector>

#include <Imagine/Features.h>

using namespace std;
using namespace Imagine;

/// Similarity transform.
class Similarity {
	FloatPoint2 center;
	float angle,zoom;
public:
	Similarity(FloatPoint2 c,float a,float z):center(c),angle(a),zoom(z) {}
	inline FloatPoint2 operator()(const FloatPoint2& P) const {
		FloatPoint2 m=(P-center)/zoom;
		double a=angle*3.14/180;
		double c=cos(a), s=sin(a);
		FloatPoint2 n=FloatPoint2(float(c*m.x()-s*m.y()),
			float(c*m.y()+s*m.x()))+center;
		return n;
	}
};

template <typename T>
Image<byte> transformImage(const Image<byte>& I,T t){
	int w=I.width(),h=I.height();
	Image<byte> O(w,h);
	for (int j=0;j<h;j++) {
		for (int i=0;i<w;i++) {
			FloatPoint2 P=t(FloatPoint2(float(i),float(j)));
			if (P.x()<0 || P.x()>w-2 || P.y()<0 || P.y()>h-2)
				O(i,j)=0;
			else
				O(i,j)=byte(I.interpolate(P));
		}
	}
	return O;
}

/// Good/bad criterion based on known transform.
template <typename T,typename S>
class TransformCrit {
const S& trans;
float thres;
public:
TransformCrit(const S& t,float th):trans(t),thres(th) {}
inline bool operator()(const Array<typename T::Feature>& F1,size_t i1,const Array<typename T::Feature>& F2,size_t i2) const {
	FloatPoint2 a=F1[i1].pos,b=F2[i2].pos;
	return norm(trans(b)-a)<thres;
}
};

/// Good/bad criterion based on fundamental matrix.
template <typename T>
class FCrit {
FMatrix<float,3,3> F;
float thres;
public:
FCrit(const FMatrix<float,3,3>& f, float th): F(f), thres(th) {}
inline bool operator()(const Array<typename T::Feature>& F1,size_t i1,const Array<typename T::Feature>& F2,size_t i2) const {
	FloatPoint2 a=F1[i1].pos,b=F2[i2].pos;
	FVector<float,3> m1(a.x(), a.y(), 1.0f), m2(b.x(), b.y(), 1.0f);
	FVector<float,3> line = F*m1;
	float d1 = abs((line*m2)/sqrt(line[0]*line[0]+line[1]*line[1]));
	line = tmult(F,m2);
	float d2 = abs((line*m1)/sqrt(line[0]*line[0]+line[1]*line[1]));
	return max(d1,d2)<thres;
}
};

template <typename T,typename Crit>
void CheckMatches(const Array<T>& F1,const Array<T>& F2,
			  const MatchList& matches,const Crit& crit) {
				  int n=0,o=0;
				  for (MatchList::const_iterator it=matches.begin();it!=matches.end();++it) {
					  if (crit(F1,it->first,F2,it->second))
						  n++;
					  else
						  o++;
				  }
				  cout << n << " correct matches and " << o << " outliers." << endl;
}

template <typename Det,typename Crit> 
void testDetector(string msg,const Det& D,
			  const Image<byte>&I1,const Image<byte>&I2,
			  const Crit& crit)
{
cout << "============ " << msg << " ================" << endl;
clearWindow();
display(I1);
display(I2,0,I1.height());
typedef typename Det::Feature Feat;
Timer tm;
tm.reset();
Array<Feat> feats1=D.run(I1);
cout << feats1.size() << " features in " << tm.lap() << " sec" << endl;
writeFeaturePoints(feats1,"keys.txt");	// Check a W/R
readFeaturePoints(feats1,"keys.txt");
drawFeatures(feats1);

tm.reset();
Array<Feat> feats2=D.run(I2);
cout << feats2.size() << " features in " << tm.lap() << " sec" << endl;
drawFeatures(feats2,IntPoint2(0,I1.height()));
MatchList matches=loweMatch(feats1,feats2,0.75,false);
CheckMatches(feats1,feats2,matches,crit);

cout << "click to see matches." << endl;
click();
display(I1);
display(I2,0,I1.height());
drawGoodMatches(feats1,feats2,matches,crit,
	IntPoint2(0,0),IntPoint2(0,I1.height()),1.0f,true);
cout << "click to finish " << msg << endl;
click();
}

void all_feat() {
// Testing invariance to similarity
cout << "TESTING WITH KNOWN TRANSFO" << endl;
cout << "==========================" << endl;
Image<byte> I1, I2;
load(I1,srcPath("toys.0000.png"));
int w=I1.width(),h=I1.height(); 
openWindow(w,2*h);
Similarity S(FloatPoint2(w/2.f,h/2.f),20,1.3f);
I2=transformImage(I1,S);
testDetector("VLFeat SIFT",SIFTDetector(),I1, I2,TransformCrit<SIFTDetector,Similarity>(S, 7.0f));
testDetector("Lowe SIFT",LoweSIFTDetector(),I1, I2, TransformCrit<LoweSIFTDetector,Similarity>(S, 7.0f));
testDetector("MSER",MSERDetector(),I1, I2, TransformCrit<MSERDetector,Similarity>(S, 7.0f));

// Testing two images with known fundamental matrix.
cout << "TESTING WITH KNOWN F matrix" << endl;
cout << "==========================" << endl;
load(I1,srcPath("toys.0000.png"));
	load(I2,srcPath("toys.0200.png"));
	ifstream str(srcPath("toys_F_0200_0000.txt"));
	FMatrix<float,3,3> F;
	str >> F;
	F=transpose(F); // From 0 to 200
	testDetector("VLFeat SIFT",SIFTDetector(),I1, I2, FCrit<SIFTDetector>(F, 10.0f));
	testDetector("Lowe SIFT",LoweSIFTDetector(),I1, I2, FCrit<LoweSIFTDetector>(F, 10.0f));
	testDetector("MSER",MSERDetector(),I1, I2, FCrit<MSERDetector>(F, 10.0f));
}

void SiftExtract(int argc,char* argv[])
{
	if (argc!=3 && argc!=4) {
		cerr << argv[0] << " image points" << endl;
		cerr << argv[0] << " image1 matches image2" << endl;
		exit(1);
	}
	Image<byte> I1,I2;
	string output(argv[2]);
	cout << "Loading " << argv[1] << endl;
	load(I1,argv[1]);
	SIFTDetector D;
	D.setFirstOctave(0); // -1: double image
	D.setNumOctaves(-1); // -1 = max
	D.setNumScales(3); // scales / octave
	cout << "Computing Sifts" << endl;
	Array<SIFTDetector::Feature> F1=D.run(I1);
	cout << "Found " << F1.size() << " points" << endl;
	if(argc == 4) {
		cout << "Loading " << argv[3] << endl;
		load(I2,argv[3]);
		Array<SIFTDetector::Feature> F2=D.run(I2);
		cout << "Found " << F2.size() << " points" << endl;
		MatchList matches=loweMatch(F1,F2,0.75,true);
		cout << matches.size() << " matches" << endl;

		cout << "Writing result to " << output << endl;
		ofstream str(output.c_str());
		MatchList::const_iterator it = matches.begin();
		for(; it != matches.end(); ++it)
			str << F1[it->first].x()  << " " << F1[it->first].y() << " "
			<< F2[it->second].x() << " " << F2[it->second].y() << endl;
	} else {
		cout << "Writing result to " << output << endl;
		writeFeaturePoints(F1,output,true);
	}
}

int main(int argc,char*argv[]) {
	if(argc > 1) // SIFT extraction
		SiftExtract(argc,argv);
	else {		// Testing all Feats
		all_feat();
		endGraphics();
	}
	return 0;
}
