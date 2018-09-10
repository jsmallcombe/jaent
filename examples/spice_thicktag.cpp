
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
// 	experiment.set_targ(target(20,40,0.6,180.0,0,40,79,197));
	experiment.set_targ(target(20,40,0.6,180.0,0,10,47,108));			

	// Set the beam
	experiment.set_beam(18,36,120);
	
	//we want an excited beam like so have to force it
// 	experiment.set_reco(34,70);
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
	experiment.set_gamma(1000.0);
// 	experiment.set_lifetime_ns(0.001,1.54,19.32);//includes densities
	experiment.set_lifetime_ns(1.001,1.54,10.5);//includes densities
	
	cout<<endl;
	experiment.print_target();
	experiment.print_decay();
	
	// Draw the target and decay things
	// quite important to check things are set and working right
// 	can_view->cd(1);
// 	gPad->cd(3);
// 	experiment.draw_decay_Z(1,false);
	
	can_view->cd();
	experiment.draw_target_interaction(0,false,500000);
// 	experiment.draw_target_interaction();


	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	