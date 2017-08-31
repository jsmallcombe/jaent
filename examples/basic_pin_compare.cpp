
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


int main(int argc, char *argv[]){	
// int argcb;
// char **argvb;	
TApplication *app = new TApplication("app", &argc, argv);
TH1::AddDirectory(kFALSE);//avoid name overwrites but could cause other memory histogram issues
//        h->SetDirectory(0);          for the current histogram h
	
	gStyle->SetOptStat(0);//Turn off annoying stats box
	exp_core experiment(60);//Create the simulation class, draw size 60 mm
	
	cout<<endl<<"START"<<endl;	
	
	TFile* output_file= new TFile("pin_compare.root","RECREATE");//Create output file
	gROOT->cd();//cd out of output file

	TCanvas * can_view;

	bool quick=true;//Does worse simulation but quick
	bool simulate=true;//Set false for just drawing when writing new geometry
	bool compare=true;//Does the S3
	
///////////////////////////////////////////////
/////////////////  S3   //////////////////////
//////////////////////////////////////////////
if(compare){
	spice_auto_setup(experiment,0,1);

	can_view = new TCanvas("S3_view", "S3_view", 1000, 600);
	can_view->Divide(2,1);
	can_view->cd(1);
	experiment.draw_exp();
	can_view->cd(2);
	experiment.draw_exp(3);
	output_file->cd();
		can_view->Write();
	gROOT->cd();
	
	if(simulate){
		TH1D h1;
		if(quick)h1= experiment.theta_cover_auto_norm(true,50000);
		else h1= experiment.theta_cover_auto_norm(true);
		h1.SetName("S3_theta");
		h1.SetTitle("S3_theta");
		h1.SetLineColor(2);
		output_file->cd();
			h1.Write();
		gROOT->cd();
	}
}

//////////////////////////////////////////////////////
/////////////////  Quadrants Parameters   /////////////
///////////////////////////////////////////////////////

TRotation rota,rotb;

double safety_gap=0.5;

double horizontal_spacing=10.6+safety_gap;
double vertical_spacing=12.4+safety_gap;
double vertical_start_offset=12.4+safety_gap;

double centre_pixel_shift=5.5;
//positive for bigger gap. Negative for overlap
double horizontal_shift=0;//+safety_gap; 
double vertical_shift=0;//+safety_gap; 

double ZZ=20;

//////////////////////////////////////////////////////
/////////////////  Quadrants   //////////////////////
///////////////////////////////////////////////////////

spice_auto_setup(experiment,0,0);

for(int k=0;k<4;k++){//quad
	for(int i=0;i<=3;i++){ //A axis (x)
		int x=abs(i);
		double XX=(horizontal_spacing+horizontal_shift)*x;
			for(int j=0;j<3;j++){ // B axis (y)
			  
		    //  cout<<endl<<k<<" "<<i<<" "<<j<<flush;
				int y=abs(j);
				double YY=vertical_start_offset+(vertical_spacing+vertical_shift)*j;
				if(i==0)YY+=centre_pixel_shift;
				if(x==3 && y==2)continue;

				TVector3 vec(-XX,YY,ZZ);
				vec*=rota;
// 				add_PIN(experiment,vec,rotb,false);
				add_PIN(experiment,vec,rotb,true);
				
			}
	}

	//rotate for next quadrant
	rota.RotateZ(pi/2);
	rotb.RotateZ(pi/2);		
}

	////// write /////
	
	can_view = new TCanvas("Quad_view", "Quad_view", 1000, 600);
	can_view->Divide(2,1);
	can_view->cd(1);
	experiment.draw_exp();
	can_view->cd(2);
	experiment.draw_exp(3);
	output_file->cd();
		can_view->Write();
	gROOT->cd();
	
//////////////////////////////////////////////////////
/////////////////  Simulation   //////////////////////
///////////////////////////////////////////////////////
	
	if(simulate){
		TH1D h1;
		if(quick)h1=experiment.theta_cover_auto_norm(true,50000);
		else h1=experiment.theta_cover_auto_norm(true);
		h1.SetName("Quad_theta");
		h1.SetTitle("Quad_theta");
		h1.SetLineColor(1);
		output_file->cd();
			h1.Write();
		gROOT->cd();
	}

		
	output_file->Write();
	output_file->Close();

	cout<<endl<<"///////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}
