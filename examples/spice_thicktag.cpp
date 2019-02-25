
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

	exp_core experiment(150);
	
	cout<<endl<<"START"<<endl;	
	
	//
	// Everything above here is basic setup
	//
	
	//
	// We are interested in a full target interaction, so we first create a "real" target
	// This target is 0.5 mg calcium and downstream 40 mg of gold backing
	//
 	experiment.set_targ(target(20,40,0.6,180.0,0,40,79,197));
// 	experiment.set_targ(target(20,40,0.6,180.0,0,10,47,108));			

	// Set the beam
	experiment.set_beam(18,36,120);
	
	//we want an excited beam like so have to force it
//  	experiment.set_reco(34,72);
	experiment.set_reco(38,76);
	
	// Recoil will be excited 
	experiment.set_E_star(1.5);
	
	// Now the primary reaction is set, print it to make sure we didnt F up.
	experiment.print_reaction();

	
	//
	// Now we set the target interaction and decay information
	//
	// setup the required thick target interaction goes roughly as E^2
	experiment.set_target_interaction(2);
	//experiment.set_gamma(1000.0);
    experiment.set_ICE(700);
 	experiment.set_lifetime_ns(0.001,1.54,19.32);//includes densities
// 	experiment.set_lifetime_ns(0.001,1.54,10.5);//includes densities
	
//     experiment.set_stopper_seperation(0.001);
    //
    //  General Checking Functions
    //
// 	cout<<endl;
// 	experiment.print_target();
// 	experiment.print_decay();
//     experiment.draw_exp();
//    experiment.draw_hits_3D(0,2,true,0.75,1,true,20);
//     experiment.draw_hit_pattern_2D(1,20000,false,0,0,0,1);
//     experiment.draw_hits_3D(1,1,false,0.75,10,true,120);

    //
    //  For looking at the kinematics of electrons
    //  
    add_S3(experiment,-12,true);
    for(int i=0;i<24;i++)experiment.set_valid_part(i,0,0,0,1);
    experiment.print_detectables();
	experiment.basic_hit_count(2000000,false,0,0,0,1);
    TH2D dethist("dethist","Electrons;Channel;Energy [keV]",24,0,24,1000,0,1000);
    for(int i=0;i<24;i++){
        detector a=experiment.get_det(i);
        TH1* h=&a.energyG;
        h=add_resolution(h,5,0.7);
        for(int j=1;j<=1000;j++){
            double c=h->GetBinContent(j*2)+h->GetBinContent(j*2-1);
            dethist.SetBinContent(i,j,c);
        }
    }
    dethist.DrawCopy("colz");
        
    
    
	// Draw the target and decay things
	// quite important to check things are set and working right
// 	can_view->cd(1);
// 	gPad->cd(3);
    
 // This is for inflight   
//  	experiment.draw_decay_Z(1,false);
	
    //This is for in target
// 	experiment.draw_target_interaction(0,false,5000000);


	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	
