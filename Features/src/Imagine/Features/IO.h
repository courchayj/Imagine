// ===========================================================================
// Imagine++ Libraries
// Copyright (C) Imagine
// For detailed information: http://imagine.enpc.fr/software
// ===========================================================================

namespace Imagine {
	/// \addtogroup Features
	/// @{

	/// Draw on FP.
	/// Draws a feature point.
	/// \tparam T FP type
	/// \param fp FP to draw
	/// \param off display offset
	/// \param c color
	/// \param drawAngle true=draw a radius showing direction
	/// \param fillIt true=draw a filled circle
	/// \param z scaling factor (of the whole image)
	/// 
	/// \dontinclude Features/test/test.cpp \skip matching()
	/// \skipline draw one feature
	template <typename T>
	void drawFeature(
		const T& fp,
		IntPoint2 off=IntPoint2(0,0),
		Color c=RED,
		bool drawAngle=true,
		bool fillIt=false,
		float z=1.f) {
			if (fillIt)
				fillCircle(IntPoint2(fp.pos*z)+off,int(fp.scale*z+1),c);
			else
				drawCircle(IntPoint2(fp.pos*z)+off,int(fp.scale*z+1),c);
			if (drawAngle)
				drawLine(IntPoint2(fp.pos*z)+off,
				IntPoint2(fp.pos*z+FVector<float,2>(float(fp.scale*cos(fp.angle)*z),float(fp.scale*sin(fp.angle)*z)))+off,
				c,1);
	}
	/// Draw FPs.
	/// Draws features points.
	/// \tparam T FP type
	/// \param fps array of FPs
	/// \param off display offset
	/// \param c color
	/// \param drawAngle true=draw a radius showing direction
	/// \param z scaling factor (of the whole image)
	/// 
	/// \dontinclude Features/test/test.cpp \skip matching()
	/// \skipline draw features
	/// \until ...
	template <typename T>
	void drawFeatures(const Array<T>& fps,
		IntPoint2 off=IntPoint2(0,0),
		Color c=RED,
		bool drawAngle=true,
		float z=1.f) {

			noRefreshPush();
			for (size_t i=0;i<fps.size();i++)
				drawFeature(fps[i],off,c,drawAngle,false,z);
			noRefreshPop();
	}

	/// Draw Matches.
	/// Draws matches between two images, one different (random) color each. Possibility to link matches with a line.
	/// \tparam T FP type
	/// \param F1,F2 arrays of FPs
	/// \param L list of matches
	/// \param off1,off2 display offsets
	/// \param z scaling factor (of the whole image)
	/// \param drawLines true=link matches with lines
	/// 
	/// \dontinclude Features/test/test.cpp \skip matching()
	/// \skipline draw features
	/// \until ...
	template <typename T>
	void drawMatches(const Array<T>& F1,const Array<T>& F2,
                     const MatchList& L,
                     IntPoint2 off1, IntPoint2 off2,
                     float z=1.f, bool drawLines=false) {
			for (MatchList::const_iterator it=L.begin();it!=L.end();++it) {
				Color c(intRandom(0,255),intRandom(0,255),intRandom(0,255));
				size_t i1=it->first;
				size_t i2=it->second;
				drawFeature(F1[i1],off1,c,false,true,z);
				drawFeature(F2[i2],off2,c,false,true,z);
                if(drawLines)
                    drawLine(IntPoint2(F1[i1].pos*z)+off1,
                             IntPoint2(F2[i2].pos*z)+off2,c,1);
			}
	}
	/// Draw Good Matches.
	/// Draws matches between two images, good ones in GREEN, bad ones in RED. good/bad refers to some given criterion!
	/// \tparam T FP type
	/// \tparam C criterion type
	/// \param F1,F2 arrays of FPs
	/// \param L list of matches
	/// \param crit good/bad criterion. \p crit(F1,i1,F2,i2)  should be true is OK, false is not
	/// \param off1,off2 display offsets
	/// \param z scaling factor (of the whole image)
	/// \param drawLines true=link matches with lines
	/// 
	/// \dontinclude Features/test/test.cpp 
	/// \skipline Good/bad criterion based on fundamental matrix
	/// \until ...
	/// \skipline draw good matches
	template <typename T,typename C>
	void drawGoodMatches(const Array<T>& F1,const Array<T>& F2,
                         const MatchList& L, C crit,
                         IntPoint2 off1, IntPoint2 off2,
                         float z=1.f, bool drawLines=false) {
			for (MatchList::const_iterator it=L.begin();it!=L.end();++it) {
				size_t i1=it->first;
				size_t i2=it->second;
				Color c=crit(F1,i1,F2,i2)?GREEN:RED;
				drawFeature(F1[i1],off1,c,false,true,z);
				drawFeature(F2[i2],off2,c,false,true,z);
                if(drawLines)
                    drawLine(IntPoint2(F1[i1].pos*z)+off1,
                             IntPoint2(F2[i2].pos*z)+off2,c,1);
			}
	}
	/// Read FPs.
	/// Reads feature points to a file.
	/// \tparam T FP type
	/// \param feats array of FPs
	/// \param name file name
	/// \param loweCompatibility true=Lowe's format (x/y swapped and rounded values)
	/// 
	/// \dontinclude Features/test/test.cpp \skip matching()
	/// \skipline read FP
	/// \until ...
	template <typename T>
	bool readFeaturePoints(Array<T>& feats,std::string name,bool loweCompatibility=false) {
		std::ifstream f(name.c_str());
		if (!f.is_open()) {
			std::cerr << "Cant open file" << std::endl;		
			return false;
		}
		int nb,sz;
		f >> nb >> sz;
		typename T::Desc foo;
		if (sz!=foo.size()) {
			std::cerr << "Bad desc size" << std::endl;		
			return false;
		}
		feats.setSize(nb);
		for (int i=0;i<nb;i++) {
			if (loweCompatibility)
				f >> feats[i].pos[1] >> feats[i].pos[0] >> feats[i].scale >> feats[i].angle;
			else
				f >> feats[i].pos[0] >> feats[i].pos[1] >> feats[i].scale >> feats[i].angle;
			for (int k=0;k<sz;k++) { 
				double x;
				f >> x;
				feats[i].desc[k]=typename T::Desc::value_type(x);
			}
		}
		return true;
	}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
	inline float make_round(float x, int n)
	{
		float p=pow(10.f,n);
		x*=(p);
		x=(x>=0)?floor(x):ceil(x);
		return x/p;			
	}
#endif

	/// Write FPs.
	/// Writes feature points to a file.
	/// \tparam T FP type
	/// \param feats array of FPs
	/// \param name file name
	/// \param loweCompatibility true=Lowe's format (x/y swapped and rounded values)
	/// 
	/// \dontinclude Features/test/test.cpp \skip matching()
	/// \skipline write FP
	/// \until ...
	template <typename T>
	bool writeFeaturePoints(Array<T>& feats,std::string name,bool loweCompatibility=false) {
		std::ofstream f(name.c_str());
		if (!f.is_open()) {
			std::cerr << "Cant open file" << std::endl;		
			return false;
		}
		typename T::Desc foo;
		int sz=foo.size();
		f << feats.size() << " " << sz << std::endl;
		for (size_t i=0;i<feats.size();i++) {
			if (loweCompatibility)
				f <<make_round(feats[i].y(),2) << ' ' << make_round(feats[i].x(),2) << ' ' << make_round(feats[i].scale,2) << ' ' << make_round(feats[i].angle,3) << std::endl;
			else 
				f << feats[i].pos[0] << ' ' << feats[i].pos[1] << ' ' << feats[i].scale << ' ' << feats[i].angle << std::endl;
			for (int k=0;k<sz;k++) {
				f << double(feats[i].desc[k]) << ' '; 
				if ((k+1)%20==0 || k==sz-1)
					f << std::endl;
			}
		}
		return true;
	}
	///@}

}
