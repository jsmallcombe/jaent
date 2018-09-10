
//
//
// James's not Geant4 : jaent V1.3.1
// james.smallcombe@outlook.com 01/4/2015
//
//


#ifndef __DETECTORCLASS_H_INCLUDED__   // if x.h hasn't been included yet...
#define __DETECTORCLASS_H_INCLUDED__   //   #define this so the compiler knows it has been included

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <vector>

// #include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRotation.h>
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TRandom2.h>
#include <TPolyLine3D.h>
// #include <range.h>
#include <TGraph.h>

// #include "james_nuclear_data_ob.h"
// #include "james_detector_func.h"
#include "james_physics.h"
// #include "james_fission.h"
// #include "james_target.h"



using namespace std;



/////////////////////////////////////////////////////////////////////////////////////
/////////// "exp_core" The core class for holding experimental setups     ///////////
///////////           and performinge experiment planning simms           ///////////  
/////////////////////////////////////////////////////////////////////////////////////
class detector
{
    public:
	    
	//////////////////////////////////////////////
	/////////// Public Data Members    ///////////
	//////////////////////////////////////////////
	    
	TH1D energy,energyG,mass;
	TH2D hit_pat,energtTheta;
	TPolyLine3D* ddraw3;
	TGraph*	  ddraw2;
	string name;
		
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	/////////// Public Members Functions  ///////////
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	
	// constructor
	detector();
	detector(const detector&);//copy constructor
	detector(vector<double>,vector<double>,TVector3=TVector3(0.0,0.0,1.0),TRotation=TRotation());
	
	//assignment operator
	detector& operator=(const detector&);
	
	// destructor
	~detector();
	
	
	TVector3 hit_pos_vec(TVector3,TVector3* =NULL);
	TVector3 hit_pos_vec(TLorentzVector*,TVector3* =NULL);
	TVector3 hit_pos_vec(double, double);

	bool hit_offset_check(TVector3,TVector3*);
	bool hit_offset_check(TLorentzVector*,TVector3*);	
	
	void set_draw3();
	void add_draw(TH3*);
	void add_draw(TH2*,int=3);
	void add_draw(int=3);

	
	TGraph* tp(){return detector_tp;};
	
	double X(){return detector_world_place.X();};
	double Y(){return detector_world_place.Y();};
	double Z(){return detector_world_place.Z();};
	
	double Xn(){return detector_norm.X();};
	double Yn(){return detector_norm.Y();};
	double Zn(){return detector_norm.Z();};
	
	void Fill(TLorentzVector*,TVector3* =NULL);
	
    private:
	    
	//////////////////////////////////////////////
	/////////// Private Data Members    //////////
	////////////////////////////////////////////// 

	//detectors
	TGraph* detector_tp;
	TGraph* detector_xy;
	TVector3 detector_world_place;
	TRotation detector_world_rotate;
	TVector3 detector_norm;	   
	double plane;

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	/////////// PRIVATE Members Functions  //////////
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

};


TH1D* add_resolution(TH1D*,double=0.46,double=122);

#endif 
