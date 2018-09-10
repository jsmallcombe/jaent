
//
//
// James's not Geant4 : jaent V1.3.1
// james.smallcombe@outlook.com 01/4/2015
//
//

#include <TApplication.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TProfile2D.h>
#include <TPaveText.h>
#include <TMacro.h>
#include <TGraph.h>
#include <TBrowser.h>

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>


#include "exp_core.h"
#include "detector_class.h"
#include "auto_setup.h"
using namespace std;

//Unused historic data from astar
TGraph P_range();
TGraph P_energy();

int main(int argc, char *argv[]){	
// int argcb;
// char **argvb;	
TApplication *app = new TApplication("app", &argc, argv);
TH1::AddDirectory(kFALSE);//avoid name overwrites but could cause other memory histogram issues
//        h->SetDirectory(0);          for the current histogram h


	gStyle->SetOptStat(0);
	exp_core expr(80.0);
	
	vector<vector<int>> BadPix={{5,0},{5,1},{5,2},{5,3},{5,4},{5,20},{5,21},{5,22},{9,23},{10,8},{15,26},{0,23},{1,23},{2,23},{3,23},{4,23},{5,23},{0,24},{1,24},{2,24},{3,24},{4,24},{5,24},{0,25},{1,25},{2,25},{3,25},{4,25},{5,25},{6,9},{6,10},{6,11},{6,12},{6,13},{6,14},{6,15},{5,5},{5,6},{4,0},{4,1},{5,19},{5,26},{5,27},{5,28},{5,29},{5,30},{5,31},{4,26},{4,27},{4,28},{4,29},{4,30},{4,31},{4,2}};	
	

	TRotation rot;
	rot.RotateZ(21.5*pi/180.);
	
// 	AddGapS3Rings(expr,TVector3(0.70,1.0,32.0),rot,0,23,19,21);
// 	expr.BuildOPINTI(20,60,22.2);
	
// 	AddS3BadPix(expr,TVector3(0.70,1.0,32.0),rot,0,5,BadPix,-4);
// 	expr.BuildOPINTI(20,60);
	
	AddGapS3Rings(expr,TVector3(-0.2,1.52,20.0),rot,0,23,5,5);
	expr.BuildOPINTI(20,60);

	


//  add_square(expr,TVector3(20,0,30),TRotation(),25.);
// 	new TCanvas("can_view", "can_view", 800, 800);
// 	expr.set_E_beam(5.0);
// 	expr.set_target_primary_offset(-4,0,-15);
// 	expr.draw_hits_3D(0,1,false,0.8,1,0);
	
	TCanvas * can_view = new TCanvas("can_view", "can_view", 1200, 600);
	
	can_view->Divide(2,1);
	can_view->cd(1);
	expr.draw_phi();
	can_view->cd(2);
	expr.draw_exp();	
	expr.ReadDrawMesh();
	
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
