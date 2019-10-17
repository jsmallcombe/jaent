
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
	
	TCanvas * can_view = new TCanvas("can_view", "can_view", 1200, 600);
	gStyle->SetOptStat(0);
	can_view->cd();
	can_view->Divide(2,2);
	exp_core experiment(150);
	
	cout<<endl<<"START"<<endl;	
	
	//
	// Everything above here is basic setup
	//    
    experiment.set_targ(target(20,40,0.01,180.0,0,40,79,197));
	experiment.set_beam(18,36,120);
	experiment.set_reco(35,73);
	experiment.set_E_star(4.0);
	experiment.print_reaction();
    
	experiment.set_target_interaction(2);
	experiment.set_gamma(1166);
 	experiment.set_lifetime_ns(0.01,1.54,19.32);//includes densities
	
	cout<<endl;
	experiment.print_target();
	
	//This line changes this to a plunger experiment
// 	experiment.set_stopper_seperation(0.01);
	
	//
	// Now the reaction is set we add out detectors
	//
	

	double theta=pi/4;
    
	for(int i=0;i<4;i++){
        double phi=pi*i/2;
		TRotation rot;
		rot.RotateY(theta);
		rot.RotateZ(phi);
		TVector3 vec(120*sin(theta)*cos(phi),120*sin(theta)*sin(phi),120*cos(theta));
		add_square(experiment,vec,rot,80);
		experiment.set_valid_part(experiment.detN()-1,0,0,0,1);
	}
	
	theta=pi/2;
	for(int i=0;i<8;i++){
        double phi=pi*i/4;
		TRotation rot;
		rot.RotateY(theta);
		rot.RotateZ(phi);
		TVector3 vec(120*sin(theta)*cos(phi),120*sin(theta)*sin(phi),120*cos(theta));
		add_square(experiment,vec,rot,80);
		experiment.set_valid_part(experiment.detN()-1,0,0,0,1);
	}
	
    add_annulus(experiment,-20,35,40);
    experiment.set_valid_part(12,0,0,0,0);
	
	experiment.print_detectables();
	experiment.print_decay();

	//Draw the geometry in several ways
	can_view->cd(1);
	experiment.draw_exp();
	can_view->cd(2);
	experiment.draw_target_interaction(0,false,50000);
// 	experiment.draw_decay_Z(1,false);

	// Next we simulate some events which will be recorded as good and saved to detector
	experiment.basic_hit_count(100000,false,0,0,0,1);
	
	// Finally draw detector energy spectra
	// We add some intrinsic resolution to the gamma data
	detector a=experiment.get_det(0);
    TH1* one=add_resolution(&a.energyG);
    
	for(int i=1;i<12;i++){
        TH1* Yi =(TH1*)experiment.get_det(i).energyG.Clone();
        a.energyG.Add(Yi);
        delete Yi;
    }
    
    TH1* sum=add_resolution(&a.energyG);
    
    can_view->cd(3);
    sum->Draw();
    sum->GetXaxis()->SetRangeUser(650,750);
    one->Draw("same");
    
    experiment.set_valid_part(12,0,0,0,1);
    
    for(int i=0;i<12;i++){
		experiment.set_valid_part(i,0,0,0,0);
    }

    
    experiment.set_ICE(700);
	experiment.basic_hit_count(1000000,false,0,0,0,1);
    
	can_view->cd(4);
	detector b=experiment.get_det(12);
    TH1* sili=add_resolution(&b.energyG);
    sili->Draw();
    sili->GetXaxis()->SetRangeUser(600,750);
    
    
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	
