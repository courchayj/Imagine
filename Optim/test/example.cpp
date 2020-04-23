#include <Imagine/Optim.h>

using namespace std;
using namespace Imagine;

// Types needed for least median squares
typedef FVector<double,2> FloatPoint2;
typedef std::pair<double,double> Line;

// Functor: estimate line equation from a sample of two points
class LineEstimator {
public:
	template <class OutputLineIterator>
	void operator() (const FArray<FloatPoint2*,2>& points, OutputLineIterator lines) const {
		Line line;
		line.first = ( points[1]->y() - points[0]->y() ) / ( points[1]->x() - points[0]->x() );
		line.second = points[0]->y() - points[0]->x() * line.first;
		*lines++ = line;
	}
};

// Functor: estimate square residual of point relatively to a line equation
class LineResidual {
public:
	double operator() (const Line& line, const FloatPoint2& point) const {
		double err = line.first * point.x() + line.second - point.y();
		return err*err;
	}
};

void testLeastMedianOfSquares() {
    cout << "-- Least median of squares --" << endl;
	initRandom(0);

	// Generate points on line y=2x+3, with Gaussian noise and 20% of outliers
	size_t size = 1000;
    cout << "Estimate noisy line y=2x+3+noise from " << size << " samples"<<endl;
	FloatPoint2* points = new FloatPoint2[size];
	for(size_t i=0;i<size;i++) {
		points[i].x() = double(i);
		if (i%5==0)
			points[i].y() = 0;
		else
			points[i].y() = 2. * double(i) + 3. + gaussianRandom()*.1;
	}

	// Robustly estimate the line equation using least median squares
	Line line;
	double outlier_thres;
	double median_res = leastMedianOfSquares<2>(points, points+size, LineEstimator(), LineResidual(), line, &outlier_thres);
	cout << "Estimated line: y = " << line.first << " * x + " << line.second << endl;
	cout << "Median square residual: " << median_res << endl;

	// Count inliers
	size_t inlier_size = 0;
	for (size_t i=0;i<size;i++)
		if ( LineResidual()(line,points[i]) < outlier_thres )
			inlier_size++;
	cout << "Number of inliers: " << inlier_size << endl;
	
	// should now indeed perfom a least suqre estimation with all inliers
	
	delete [] points;
}

void testRansac() {
    cout << "-- RANSAC --" << endl;
	initRandom(0);

	// Generate points on line y=2x+3, with Gaussian noise and 20% of outliers
	size_t size = 1000;
    cout << "Estimate noisy line y=2x+3+noise from " << size << " samples"<<endl;
	FloatPoint2* points = new FloatPoint2[size];
	for(size_t i=0;i<size;i++) {
		points[i].x() = double(i);
		if (i%5==0)
			points[i].y() = 0;
		else
			points[i].y() = 2. * double(i) + 3. + gaussianRandom()*.1;
	}

	// Robustly estimate the line equation using RANSAC
	Line line;
	double outlier_thres=.2;
	size_t innb = ransac<2>(points, points+size, LineEstimator(), LineResidual(), line, outlier_thres);
	cout << "Estimated line: y = " << line.first << " * x + " << line.second << endl;
	cout << "NB of inliers: " << innb << endl;
	// should now indeed perfom a least suqre estimation with all inliers
	delete [] points;
}

/////////////////////////////////////////////////////////////////////////////
// Example illustrating the use of maxFlow (Vlad. Kolmogorov)
// This section shows how to use the library to compute a minimum cut on the following graph:
/////////////////////////////////////////////////////////////////////////////
//
//		        SOURCE
//		       /       \
//		     1/         \6
//		     /      4    \
//		   node0 -----> node1
//		     |   <-----   |
//		     |      3     |
//		     \            /
//		     5\          /1
//		       \        /
//		          SINK
//
///////////////////////////////////////////////////

void testBinaryGraphCut()
{
	Graph<int,int,int> g(/*estimated # of nodes*/ 2, /*estimated # of edges*/ 1); 

	g.addNode(2); 

	g.addTweights( 0,   /* capacities */  1, 5 );
	g.addTweights( 1,   /* capacities */  6, 1 );
	g.addEdge( 0, 1,    /* capacities */  4, 3 );

	int flow = g.maxFlow();

	cout << "Flow = " << flow << endl;
	for (int i=0;i<2;i++)
		if (g.whatSegment(i) == Graph<int,int,int>::SOURCE)
			cout << i << " is in the SOURCE set" << endl;
		else
			cout << i << " is in the SINK set" << endl;
}

/////////////////////////////////////////////////////////////////////////////
// Example illustrating the use of GCoptimization (Olga Veksler)
/////////////////////////////////////////////////////////////////////////////
//
//  Optimization problem:
//  is a set of sites (pixels) of width 10 and hight 5. Thus number of pixels is 50
//  grid neighborhood: each pixel has its left, right, up, and bottom pixels as neighbors
//  7 labels
//  Data costs: D(pixel,label) = 0 if pixel < 25 and label = 0
//            : D(pixel,label) = 10 if pixel < 25 and label is not  0
//            : D(pixel,label) = 0 if pixel >= 25 and label = 5
//            : D(pixel,label) = 10 if pixel >= 25 and label is not  5
// Smoothness costs: V(p1,p2,l1,l2) = min( (l1-l2)*(l1-l2) , 4 )
// Below in the main program, we illustrate different ways of setting data and smoothness costs
// that our interface allow and solve this optimizaiton problem

// For most of the examples, we use no spatially varying pixel dependent terms. 
// For some examples, to demonstrate spatially varying terms we use
// V(p1,p2,l1,l2) = w_{p1,p2}*[min((l1-l2)*(l1-l2),4)], with 
// w_{p1,p2} = p1+p2 if |p1-p2| == 1 and w_{p1,p2} = p1*p2 if |p1-p2| is not 1

typedef short ValType;
typedef int TotValType;
//Possible choices:(short,int), (int,int) (float,float) and (double,double)

struct ForDataFn{
	int numLab;
	ValType *data;
};


ValType smoothFn(int, int, int l1, int l2)
{
	if ( (l1-l2)*(l1-l2) <= 4 ) return ValType((l1-l2)*(l1-l2));
	else return ValType(4);
}

ValType smoothFn2(int p1, int p2, int l1, int l2)
{
	return ValType(min((l1-l2)*(l1-l2),4)*((abs(p1-p2)==1)?p1+p2:p1*p2));
}

ValType dataFn(int p, int l, void *data)
{
	ForDataFn *myData = (ForDataFn *) data;
	int numLab = myData->numLab;
	
	return( myData->data[p*numLab+l] );
}



////////////////////////////////////////////////////////////////////////////////
// smoothness and data costs are set up one by one, individually
// grid neighborhood structure is assumed
//
void testGridGraphIndividually(int width,int height,int numPixels,int numLabels)
{

	int *result = new int[numPixels];   // stores result of optimization



	try{
		GCoptimizationGridGraph<ValType,TotValType> *gc = new GCoptimizationGridGraph<ValType,TotValType>(width,height,numLabels);

		// first set up data costs individually
		for ( int i = 0; i < numPixels; i++ )
			for (int l = 0; l < numLabels; l++ )
				if (i < 25 ){
					if(  l == 0 ) gc->setDataCost(i,l,(ValType)0);
					else gc->setDataCost(i,l,(ValType)10);
				}
				else {
					if(  l == 5 ) gc->setDataCost(i,l,(ValType)0);
					else gc->setDataCost(i,l,(ValType)10);
				}

		// next set up smoothness costs individually
		for ( int l1 = 0; l1 < numLabels; l1++ )
			for (int l2 = 0; l2 < numLabels; l2++ ){
				ValType cost = ValType((l1-l2)*(l1-l2) <= 4  ? (l1-l2)*(l1-l2):4);
				gc->setSmoothCost(l1,l2,cost); 
			}

		cout << "Before optimization energy is " << gc->computeEnergy();
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		cout << "\nAfter optimization energy is " << gc->computeEnergy();

		for ( int  i = 0; i < numPixels; i++ )
			result[i] = gc->whatLabel(i);

		delete gc;
	}
	catch (GCException e){
		e.report();
	}

	delete [] result;
}

////////////////////////////////////////////////////////////////////////////////
// in this version, set data and smoothness terms using arrays
// grid neighborhood structure is assumed
//
void testGridGraphDArraySArray(int width,int height,int numPixels,int numLabels)
{

	int *result = new int[numPixels];   // stores result of optimization

	// first set up the array for data costs
	ValType *data = new ValType[numPixels*numLabels];
	for ( int i = 0; i < numPixels; i++ )
		for (int l = 0; l < numLabels; l++ )
			if (i < 25 ){
				if(  l == 0 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
			else {
				if(  l == 5 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
	// next set up the array for smooth costs
	ValType *smooth = new ValType[numLabels*numLabels];
	for ( int l1 = 0; l1 < numLabels; l1++ )
		for (int l2 = 0; l2 < numLabels; l2++ )
			smooth[l1+l2*numLabels] = ValType((l1-l2)*(l1-l2) <= 4  ? (l1-l2)*(l1-l2):4);


	try{
		GCoptimizationGridGraph<ValType,TotValType> *gc = new GCoptimizationGridGraph<ValType,TotValType>(width,height,numLabels);
		gc->setDataCost(data);
		gc->setSmoothCost(smooth);
		cout << "Before optimization energy is " << gc->computeEnergy();
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		cout << "\nAfter optimization energy is " << gc->computeEnergy();
		
		for ( int  i = 0; i < numPixels; i++ )
			result[i] = gc->whatLabel(i);

		delete gc;
	}
	catch (GCException e){
		e.report();
	}

	delete [] result;
	delete [] smooth;
	delete [] data;

}
////////////////////////////////////////////////////////////////////////////////
// in this version, set data and smoothness terms using arrays
// grid neighborhood structure is assumed
//
void testGridGraphDfnSfn(int width,int height,int numPixels,int numLabels)
{

	int *result = new int[numPixels];   // stores result of optimization

	// first set up the array for data costs
	ValType *data = new ValType[numPixels*numLabels];
	for ( int i = 0; i < numPixels; i++ )
		for (int l = 0; l < numLabels; l++ )
			if (i < 25 ){
				if(  l == 0 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
			else {
				if(  l == 5 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}


	try{
		GCoptimizationGridGraph<ValType,TotValType> *gc = new GCoptimizationGridGraph<ValType,TotValType>(width,height,numLabels);

		// set up the needed data to pass to function for the data costs
		ForDataFn toFn;
		toFn.data = data;
		toFn.numLab = numLabels;

		gc->setDataCost(&dataFn,&toFn);

		// smoothness comes from function pointer
		gc->setSmoothCost(&smoothFn2);

		cout << "Before optimization energy is " << gc->computeEnergy();
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		cout << "\nAfter optimization energy is " << gc->computeEnergy();
		
		cout << endl;
		for ( int  i = 0; i < numPixels; i++ ) {
			result[i] = gc->whatLabel(i);
			cout << result[i] << " ";
		}

		delete gc;
	}
	catch (GCException e){
		e.report();
	}

	delete [] result;
	delete [] data;

}
////////////////////////////////////////////////////////////////////////////////
// Uses spatially varying smoothness terms. VH pattern
void testGridGraphDArraySArraySpatVH(int width,int height,int numPixels,int numLabels)
{
	int *result = new int[numPixels];   // stores result of optimization

	// first set up the array for data costs
	ValType *data = new ValType[numPixels*numLabels];
	for ( int i = 0; i < numPixels; i++ )
		for (int l = 0; l < numLabels; l++ )
			if (i < 25 ){
				if(  l == 0 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
			else {
				if(  l == 5 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
	// next set up the array for smooth costs
	ValType *smooth = new ValType[numLabels*numLabels];
	for ( int l1 = 0; l1 < numLabels; l1++ )
		for (int l2 = 0; l2 < numLabels; l2++ )
			smooth[l1+l2*numLabels] = ValType((l1-l2)*(l1-l2) <= 4  ? (l1-l2)*(l1-l2):4);

	// next set up spatially varying arrays V and H

	ValType *V = new ValType[numPixels];
	ValType *H = new ValType[numPixels];

	
	for ( int i = 0; i < numPixels; i++ ){
		H[i] = ValType(i+(i+1)%3);
		V[i] = ValType(i*(i+width)%7);
	}


	try{
		GCoptimizationGridGraph<ValType,TotValType> *gc = new GCoptimizationGridGraph<ValType,TotValType>(width,height,numLabels);
		gc->setDataCost(data);
		gc->setSmoothCostVH(smooth,V,H);
		cout << "Before optimization energy is " << gc->computeEnergy();
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		cout << "\nAfter optimization energy is " << gc->computeEnergy();
		
		for ( int  i = 0; i < numPixels; i++ )
			result[i] = gc->whatLabel(i);

		delete gc;
	}
	catch (GCException e){
		e.report();
	}

	delete [] result;
	delete [] smooth;
	delete [] data;


}

////////////////////////////////////////////////////////////////////////////////
// in this version, set data and smoothness terms using arrays
// grid neighborhood is set up "manually"
//
void testGeneralGraphDArraySArray(int width,int height,int numPixels,int numLabels)
{

	int *result = new int[numPixels];   // stores result of optimization

	// first set up the array for data costs
	ValType *data = new ValType[numPixels*numLabels];
	for ( int i = 0; i < numPixels; i++ )
		for (int l = 0; l < numLabels; l++ )
			if (i < 25 ){
				if(  l == 0 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
			else {
				if(  l == 5 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
	// next set up the array for smooth costs
	ValType *smooth = new ValType[numLabels*numLabels];
	for ( int l1 = 0; l1 < numLabels; l1++ )
		for (int l2 = 0; l2 < numLabels; l2++ )
			smooth[l1+l2*numLabels] = ValType((l1-l2)*(l1-l2) <= 4  ? (l1-l2)*(l1-l2):4);


	try{
		GCoptimizationGeneralGraph<ValType,TotValType> *gc = new GCoptimizationGeneralGraph<ValType,TotValType>(numPixels,numLabels);
		gc->setDataCost(data);
		gc->setSmoothCost(smooth);

		// now set up a grid neighborhood system
		// first set up horizontal neighbors
		for (int y = 0; y < height; y++ )
			for (int  x = 1; x < width; x++ )
				gc->setNeighbors(x+y*width,x-1+y*width);

		// next set up vertical neighbors
		for (int y = 1; y < height; y++ )
			for (int  x = 0; x < width; x++ )
				gc->setNeighbors(x+y*width,x+(y-1)*width);

		cout << "Before optimization energy is " << gc->computeEnergy();
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		cout << "\nAfter optimization energy is " << gc->computeEnergy();
		
		for ( int  i = 0; i < numPixels; i++ )
			result[i] = gc->whatLabel(i);

		delete gc;
	}
	catch (GCException e){
		e.report();
	}

	delete [] result;
	delete [] smooth;
	delete [] data;

}

////////////////////////////////////////////////////////////////////////////////
// in this version, set data and smoothness terms using arrays
// grid neighborhood is set up "manually". Uses spatially varying terms. Namely
// V(p1,p2,l1,l2) = w_{p1,p2}*[min((l1-l2)*(l1-l2),4)], with 
// w_{p1,p2} = p1+p2 if |p1-p2| == 1 and w_{p1,p2} = p1*p2 if |p1-p2| is not 1

void testGeneralGraphDArraySArraySpatVarying(int width,int height,int numPixels,int numLabels)
{
	int *result = new int[numPixels];   // stores result of optimization

	// first set up the array for data costs
	ValType *data = new ValType[numPixels*numLabels];
	for ( int i = 0; i < numPixels; i++ )
		for (int l = 0; l < numLabels; l++ )
			if (i < 25 ){
				if(  l == 0 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
			else {
				if(  l == 5 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
	// next set up the array for smooth costs
	ValType *smooth = new ValType[numLabels*numLabels];
	for ( int l1 = 0; l1 < numLabels; l1++ )
		for (int l2 = 0; l2 < numLabels; l2++ )
			smooth[l1+l2*numLabels] = ValType((l1-l2)*(l1-l2) <= 4  ? (l1-l2)*(l1-l2):4);


	try{
		GCoptimizationGeneralGraph<ValType,TotValType> *gc = new GCoptimizationGeneralGraph<ValType,TotValType>(numPixels,numLabels);
		gc->setDataCost(data);
		gc->setSmoothCost(smooth);

		// now set up a grid neighborhood system
		// first set up horizontal neighbors
		for (int y = 0; y < height; y++ )
			for (int  x = 1; x < width; x++ ){
				int p1 = (x-1+y*width);
				int p2 = (x+y*width);
				gc->setNeighbors(p1,p2,ValType(p1+p2));
			}

		// next set up vertical neighbors
		for (int y = 1; y < height; y++ )
			for (int  x = 0; x < width; x++ ){
				int p1 = (x+(y-1)*width);
				int p2 = (x+y*width);
				gc->setNeighbors(p1,p2,ValType(p1*p2));
			}

		cout << "Before optimization energy is " << gc->computeEnergy();
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		cout << "\nAfter optimization energy is " << gc->computeEnergy();
		
		cout << endl;
		for ( int  i = 0; i < numPixels; i++ ) {
			result[i] = gc->whatLabel(i);
			cout << result[i] << " ";
		}
		delete gc;
	}
	catch (GCException e){
		e.report();
	}

	delete [] result;
	delete [] smooth;
	delete [] data;


}

////////////////////////////////////////////////////////////////////////////////
// in this version, set data and smoothness terms using arrays
// grid neighborhood is set up "manually". 
// Uses an fn for spatially varying terms


void testGeneralGraphDArraySfn(int width,int height,int numPixels,int numLabels)
{
	int *result = new int[numPixels];   // stores result of optimization

	// first set up the array for data costs
	ValType *data = new ValType[numPixels*numLabels];
	for ( int i = 0; i < numPixels; i++ )
		for (int l = 0; l < numLabels; l++ )
			if (i < 25 ){
				if(  l == 0 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
			else {
				if(  l == 5 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}

	try{
		GCoptimizationGeneralGraph<ValType,TotValType> *gc = new GCoptimizationGeneralGraph<ValType,TotValType>(numPixels,numLabels);
		gc->setDataCost(data);
		gc->setSmoothCost(&smoothFn2);

		// now set up a grid neighborhood system
		// first set up horizontal neighbors
		for (int y = 0; y < height; y++ )
			for (int  x = 1; x < width; x++ ){
				int p1 = (x-1+y*width);
				int p2 = (x+y*width);
				gc->setNeighbors(p1,p2);
			}

		// next set up vertical neighbors
		for (int y = 1; y < height; y++ )
			for (int  x = 0; x < width; x++ ){
				int p1 = (x+(y-1)*width);
				int p2 = (x+y*width);
				gc->setNeighbors(p1,p2);
			}

		cout << "Before optimization energy is " << gc->computeEnergy();
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		cout << "\nAfter optimization energy is " << gc->computeEnergy();
		
		cout << endl;
		for ( int  i = 0; i < numPixels; i++ ) {
			result[i] = gc->whatLabel(i);
			cout << result[i] << " ";
		}

		delete gc;
	}
	catch (GCException e){
		e.report();
	}

	delete [] result;
	delete [] data;


}

////////////////////////////////////////////////////////////////////////////////
// in this version, set data and smoothness terms using arrays
// grid neighborhood is set up "manually". 
// Uses a functor to memorize spatially varying terms

// V(p1,p2,l1,l2) = w_{p1,p2}*[min((l1-l2)*(l1-l2),4)], with 
// w_{p1,p2} = p1+p2 if |p1-p2| == 1 and w_{p1,p2} = p1*p2 if |p1-p2| is not 1ValType smoothFn2(int p1, int p2, int l1, int l2)

class SF2:public GCoptimization<ValType,TotValType>::SmoothCostFunctor {
	MultiArray<ValType,4>  V;
public:
	SF2(int nodes,int labels) {
		FVector<int,4> c;
		c[0]=c[1]=nodes;
		c[2]=c[3]=labels;
		V=MultiArray<ValType,4>(c);
	}
	void setV(int s1, int s2, int l1, int l2,ValType v) {
		FVector<int,4> c;
		c[0]=s1;c[1]=s2;
		c[2]=l1;c[3]=l2;
		V(c)=v;
	}
	ValType compute(int s1, int s2, int l1, int l2) {
		FVector<int,4> c;
		c[0]=s1;c[1]=s2;
		c[2]=l1;c[3]=l2;
		return V(c);
	}
};


void testGeneralGraphDArraySfunctor(int width,int height,int numPixels,int numLabels)
{
	int *result = new int[numPixels];   // stores result of optimization

	// first set up the array for data costs
	ValType *data = new ValType[numPixels*numLabels];
	for ( int i = 0; i < numPixels; i++ )
		for (int l = 0; l < numLabels; l++ )
			if (i < 25 ){
				if(  l == 0 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}
			else {
				if(  l == 5 ) data[i*numLabels+l] = 0;
				else data[i*numLabels+l] = 10;
			}

	try{
		GCoptimizationGeneralGraph<ValType,TotValType> *gc = new GCoptimizationGeneralGraph<ValType,TotValType>(numPixels,numLabels);
		gc->setDataCost(data);
		SF2 sf2(numPixels,numLabels);
		
		// now set up a grid neighborhood system
		// first set up horizontal neighbors
		for (int y = 0; y < height; y++ )
			for (int  x = 1; x < width; x++ ){
				int p1 = (x-1+y*width);
				int p2 = (x+y*width);
				gc->setNeighbors(p1,p2);
				for ( int l1 = 0; l1 < numLabels; l1++ )
					for (int l2 = 0; l2 < numLabels; l2++ ) {
						sf2.setV(p1,p2,l1,l2,smoothFn2(p1,p2,l1,l2));
						sf2.setV(p2,p1,l1,l2,smoothFn2(p1,p2,l1,l2));
					}
			}

		// next set up vertical neighbors
		for (int y = 1; y < height; y++ )
			for (int  x = 0; x < width; x++ ){
				int p1 = (x+(y-1)*width);
				int p2 = (x+y*width);
				gc->setNeighbors(p1,p2);
				for ( int l1 = 0; l1 < numLabels; l1++ )
					for (int l2 = 0; l2 < numLabels; l2++ ) {
						sf2.setV(p1,p2,l1,l2,smoothFn2(p1,p2,l1,l2));
						sf2.setV(p2,p1,l1,l2,smoothFn2(p1,p2,l1,l2));
					}
			}

		gc->setSmoothCostFunctor(&sf2);

		cout << "Before optimization energy is " << gc->computeEnergy();
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		cout << "\nAfter optimization energy is " << gc->computeEnergy();
		
		cout << endl;
		for ( int  i = 0; i < numPixels; i++ ) {
			result[i] = gc->whatLabel(i);
			cout << result[i] << " ";
		}

		delete gc;
	}
	catch (GCException e){
		e.report();
	}

	delete [] result;
	delete [] data;


}

////////////////////////////////////////////////////////////////////////////////

void testAlphaExpansions()
{
	int width = 10;
	int height = 5;
	int numPixels = width*height;
	int numLabels = 7;

	// smoothness and data costs are set up one by one, individually
	cout << "Grid. Individual" << endl;
	testGridGraphIndividually(width,height,numPixels,numLabels);
	cout << endl;

	// smoothness and data costs are set up using arrays
	cout << "Grid. Array" << endl;
	testGridGraphDArraySArray(width,height,numPixels,numLabels);
	cout << endl;
	
	// smoothness and data costs are set up using arrays. 
	//Will pretend our graph is 
	//general, and set up a neighborhood system
	// which actually is a grid
	cout << "General. Array" << endl;
	testGeneralGraphDArraySArray(width,height,numPixels,numLabels);
	cout << endl;

	// spatially varying terms are present. VH pattern
	cout << endl;
	cout << "Grid. VH pattern" << endl;
	testGridGraphDArraySArraySpatVH(width,height,numPixels,numLabels);
	cout << endl;


	// smoothness and data costs are set up using functions
	cout << endl;
	cout << "Grid. Functions" << endl;
	testGridGraphDfnSfn(width,height,numPixels,numLabels);
	cout << endl;

	//Will pretend our graph is general, and set up a neighborhood system
	// which actually is a grid. Also uses spatially varying terms
	cout << "General. w(p,q)*f(lp,lq)" << endl;
	testGeneralGraphDArraySArraySpatVarying(width,height,numPixels,numLabels);
	cout << endl;

	//Will pretend our graph is general, and set up a neighborhood system
	// which actually is a grid. 
	// Also uses spatially varying terms via a function
	cout << "General. Function" << endl;
	testGeneralGraphDArraySfn(width,height,numPixels,numLabels);
	cout << endl;

	// Idem Functor
	cout << "General. Functor" << endl;
	testGeneralGraphDArraySfunctor(width,height,numPixels,numLabels);
	cout << endl;
}

int main()
{
    testRansac();	// RANSAC
    testLeastMedianOfSquares();	// Least median of squares
	testBinaryGraphCut();	// Maxflow (Vlad Kolm.)
	testAlphaExpansions();	// GC MRF optimization using alpha expansion (Olga Veksler)
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////
