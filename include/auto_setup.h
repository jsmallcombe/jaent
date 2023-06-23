
//
//
// James's not Geant4 : jaent V1.3.1
// james.smallcombe@outlook.com 01/4/2015
//
//


#ifndef __DET_AUTO_H_INCLUDED__   // if x.h hasn't been included yet...
#define __DET_AUTO_H_INCLUDED__   //   #define this so the compiler knows it has been included

#include "exp_core.h"
#include "detector_class.h"

using namespace std;


void mwpc_auto_setup(exp_core&,int=0,double=0,bool=false);

void add_rutherford_monitor(exp_core&,double=-36,double=176.363,double=1.2);
void gen_annulus(vector<double>&,vector<double>&,double,double);
void add_annulus(exp_core&,double=20,double=10,double=50);
void add_annulus(exp_core&,double,double,TVector3,TRotation=TRotation());

void add_micron_S(exp_core &,double,bool=false,int=3);
void add_S1(exp_core &,double,bool=true);
void add_S2(exp_core &,double,bool=true);
void add_S3(exp_core &,double,bool=true);
void add_PIN(exp_core &,TVector3,TRotation,bool=false,bool=true);
void add_square(exp_core &,TVector3,TRotation,double);

void add_target_ladder(exp_core &,double=15,double=100,double=5);

void add_chamber_cylinder(exp_core &,double=300,double=100,bool=false);
void add_chamber_cubeoid(exp_core &,double=500,double=50,double=200);
	
	
void spice_auto_setup(exp_core &,int=0,int=1,double=0);

void add_PDA(exp_core &,double=28,bool=false);



void AddGapS3Rings(exp_core &exp,TVector3 pos,TRotation rot,int InnerRing,int OuterRing,int OffSectorStart=-1,int OffSectorEnd=-1);

void AddS3BadPix(exp_core &exp,TVector3 pos,TRotation rot,int InnerRing,int OuterRing,vector<vector<int>> BadPix,int SectorOffsetN=0);





	
#endif 
