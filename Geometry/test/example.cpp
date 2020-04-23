
#include <Imagine/Geometry.h>
#include <sstream>

using namespace Imagine;
using namespace std;


typedef FVector<float,4> Point4;
typedef FVector<float,6> Point6;
typedef vector<Point4> ListPoints;
typedef vector<Point6> ListPoints6;

void epipolar(string n1,string c1,string n2,string c2,float disto=0,int fact=1)
{
	// Load two images
	Image<byte,2> I1,I2;
	if (!load(I1,n1) || !load(I2,n2)) return;

	// Load associated cameras
	Camera<float> camera1;
	if (!loadText(camera1,c1)) return;
	Camera<float> camera2;
	if (!loadText(camera2,c2)) return;

	// Correct radial distortion if needed
	Image<float,2> I1u,I2u;
	if (disto!=0) {
		camera1.distortionCoeffs()[0] = disto;
		camera1.distortionDirect() = false;
		camera2.distortionCoeffs()[0] = disto;
		camera2.distortionDirect() = false;
		I1u=correctRadialDistortion(I1,camera1);
		I2u=correctRadialDistortion(I2,camera2);
	} else {
		I1u=I1;
		I2u=I2;
	}
	// Epipolar geometry
	int w=I1.width(),h=I1.height();
	WindowType types[2]={WINDOW_2D,WINDOW_2D};
	string names[2]={"I1","I2"};
	Window W;
	W=openComplexWindow(w/fact,h/fact,"Epipolar lines",2,types,names);
	setActiveWindow(W,0);display(grey(reduce(I1u,fact)));
	setActiveWindow(W,1);display(grey(reduce(I2u,fact)));
	Camera<float>::Matrix3x3 F=camera2.fundamental(camera1);
	IntPoint2 p;
	Window wn;
	int sw;
	while (anyGetMouse(p,wn,sw)!=3) {
		setActiveWindow(wn,sw);
		drawCircle(p,4,RED);
		Camera<float>::Vector3 m1(float(p.x()*fact),float(p.y()*fact),1);
		Camera<float>::Vector3 m2=((sw==0)?F:transpose(F))*m1;
		setActiveWindow(W,1-sw);
		drawLine(int(-m2.z()/m2.x()/fact),0,
                 int((-m2.z()-m2.y()*h)/m2.x()/fact),h/fact,RED);
	}
	closeWindow(W);
}


void triangulate()
{
    cout << "-- Test triangulation --" << endl;
	Camera<float> camera1,camera2;
	if (!loadText(camera1,srcPath("toyscam0000.txt"))) return;
	if (!loadText(camera2,srcPath("toyscam0200.txt"))) return;

	FVector<float,2> p1(525,286), p2(467,365);
	FVector<float,3> P;

	Timer timer;

	Triangulation<float> triangulation;
	for (int i=0;i<100000;i++)
	{
		triangulation.clear();		
		triangulation.add(camera1,p1);
		triangulation.add(camera2,p2);
		P = triangulation.compute();
	}

	cout << timer.lap() << " seconds for 100000 triangulations" << endl;
	cout << "Position: " << P << ", Error: " << triangulation.error() << endl;
}


/// 7- or 8- point RANSAC estimation of fundamental matrix.
template <int NumPoints>
void show_fundamental(const ListPoints& lstPoints, int w1, int w2)
{
	// Run epipolar estimation
	Fundamental<float> F;
	double outlier_thresh=1.5f;
	double max_ratio = .4;

    FResidual<float> residual;
	ransac<NumPoints>(lstPoints.begin(), lstPoints.end(), 
                                 FEstimator<float>(),
                                 residual,
                                 F, outlier_thresh, max_ratio);

    // Show inliers
	int inlier_size = 0;
	for (ListPoints::const_iterator iter=lstPoints.begin(); iter != lstPoints.end();iter++)
		if (residual(F,*iter) < outlier_thresh) {
			Color c((byte)intRandom(0,255),
                    (byte)intRandom(0,255),
                    (byte)intRandom(0,255));
			fillCircle( (int) (*iter)[0] , (int) (*iter)[1],4,c);
			fillCircle(w1 + (int)(*iter)[2] ,(int) (*iter)[3],4,c);	
			inlier_size++;
		} else {
			drawCircle( (int) (*iter)[0] , (int) (*iter)[1],2,RED);
			drawCircle(w1 + (int) (*iter)[3] ,(int) (*iter)[2],2,RED);	
		}
    cout << NumPoints << "-point algorithm" << endl;
    cout << "Number of Epipolar inliers: " << inlier_size << endl;
    cout << "Inliers=full circle, outliers=empty circle" << endl;

    // Interactive display of epipolar lines
    cout << "click in an image to show epipolar line in other image"<<endl;
    cout << "Right-click to finish" << endl;
    IntPoint2 P;
    while(getMouse(P) != 3) {
        drawCircle(P,2,BLUE);
        if(P[0] < w1) { // click in left image
            FVector<float,3> l = F.rightLine((float)P[0], (float)P[1]);
            drawLine(w1,int(-l[2]/l[1]),
                     w2+w1,int((-l[2]-w2*l[0])/l[1]),BLUE);
        } else { // click in right image
            FVector<float,3> l = F.leftLine((float)P[0]-w1,(float)P[1]);
            drawLine(0,int(-l[2]/l[1]),
                     w1,int((-l[2]-w1*l[0])/l[1]),BLUE);
        }
    }
}

// Read matches file
void read_matches(ListPoints& lstPoints, const char* name)
{
	ifstream f(name);
	if(!f.is_open())
	{
		cout<<"Impossible d'ouvrir le fichier en lecture !"<<endl;
		return;
	}
    while( f.good() ) {
        std::string str;
        std::getline(f, str);
        if( f.good() ) {
            std::istringstream s(str);
            Point4 p;
            s >> p;
            if(!s.fail())
                lstPoints.push_back( p );
        }
    }
}

void ransac_fundamental() {
	cout << "-- ransac F --" << endl;
	Image<byte,2> I1,I2;
	ListPoints lstPoints;
    int w1,w2,h1,h2;

    lstPoints.clear();
	load(I1,srcPath("city0.jpg")); I1=reduce(I1,2);
	load(I2,srcPath("city1.jpg")); I2=reduce(I2,2);
	w1=I1.width(); w2=I2.width(); h1=I1.height(); h2=I2.height();
	Window W = openWindow(w1+w2,max(h1,h2));
	display(grey(I1));
	display(grey(I2),w1,0);

    read_matches(lstPoints, srcPath("citymatch01.txt"));
    show_fundamental<7>(lstPoints, w1, w2);
	display(grey(I1)); display(grey(I2),w1,0);
    show_fundamental<8>(lstPoints, w1, w2);
    closeWindow(W);
}

void ransac_trifocal() {
	cout << "-- ransac trifocal --" << endl;
	float x1, x2, y1, y2, x3, y3;
	int n,record_no,i;
	char c;
	int inlier_size = 0;
	ListPoints6 lstPoints6;
	Point6 point6;
	ifstream f1(srcPath("1.dat"));
	ifstream f2(srcPath("2.dat"));
	ifstream f3(srcPath("3.dat"));
	if(!f1.is_open())
	{cout<<"Impossible d'ouvrir le fichier en lecture !"<<endl;
	return;}

	record_no = 167;
	for (i=0; i<(record_no ); i++) {
		f1>>x1>>c>>y1>>c>>n>>c>>n;
		f2>>x2>>c>>y2>>c>>n>>c>>n;
		f3>>x3>>c>>y3>>c>>n>>c>>n;
		point6[0] = (float)x1;
		point6[1] = (float)y1;
		point6[2] = (float)x2;
		point6[3] = (float)y2;
		point6[4] = (float)x3;
		point6[5] = (float)y3;
		lstPoints6.push_back(point6);
	}
	f1.close();
	f2.close();
	f3.close();
	float out_t = 2.24f;
	Trifocal<float,7,ListPoints6::const_iterator> trif;

    ransac<7>(lstPoints6.begin(), lstPoints6.end(),
                         TrifEstimator<float,ListPoints6::const_iterator>(),
                         TrifResidual<float,ListPoints6::const_iterator>(),
                         trif,
                         out_t);

	inlier_size = 0;
	for (ListPoints6::const_iterator iter=lstPoints6.begin(); iter != lstPoints6.end();iter++)
		if ( TrifResidual<float,ListPoints6::const_iterator>()(trif,*iter) < out_t)
		{	
			inlier_size++;
		}
		cout << "Number of Trifocal inliers on 167 : " << inlier_size << endl;
}

int main(int argc, char** argv) {
	if(argc != 1) { // TestEpi binary
        if (argc!=5 && argc!=6) {
            cerr << "Usage: " << argv[0] << " ima1 cam1 ima2 cam2  [fact]"
                 << endl;
            exit(1);
        }
        epipolar(argv[1],argv[2],argv[3],argv[4],0,(argc==6)?atoi(argv[5]):1);
    } else {
        cout << "-- Test F from camera matrices --" << endl;
        cout << "click in image to show epipolar line in other image" << endl;
        cout << "Right-click to finish" << endl;
        epipolar(srcPath("toys.0000.png"),srcPath("toyscam0000.txt"),
                 srcPath("toys.0200.png"),srcPath("toyscam0200.txt"));
		triangulate();
		ransac_fundamental();
		ransac_trifocal();
    }
    endGraphics();
	return 0;
}
