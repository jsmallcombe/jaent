
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
	can_view->Divide(2,1);
	exp_core experiment(150);
	
	cout<<endl<<"START"<<endl;	
	
	//
	// Everything above here is basic setup
	//
	
	//
	// We are interested in a full target interaction, so we first create a "real" target
	// This target is 2.165 um carbon and downstream 14.917 um of gold backing
	//
	experiment.set_targ(target(6,12,giveme_areal(2.25,2.165),180.0,0,giveme_areal(19.32,14.917),79,197));		

	// Set the beam
	experiment.set_beam(36,86,256.71);
	
	//we want an excited beam like so have to force it
	experiment.set_reco(36,86);
	
	// Recoil will be excited 
	experiment.set_E_star(1.564);
	
	// Now the primary reaction is set, print it to make sure we didnt F up.
	experiment.print_reaction();
	
	//cut down on wasted simulations
	experiment.set_uniform(0,0.5*pi);
	//set_rutherford(double thetamin=0.1,double thetamax=pi)
	// experiment.reverse_primary_dist();//as we made our recoil beam-like	

	
	//
	// Now we set the target interaction and decay information
	//
	// setup the required thick target interaction goes roughly as E^2
	experiment.set_target_interaction(2);
	experiment.set_gamma(1564);
	//experiment.set_lifetime_ns(0.000426,2.25,19.32);//includes densities
	
	cout<<endl;
	experiment.print_target();
	
	//This line changes this to a plunger experiment
	//experiment.set_stopper_seperation(0.1);
	
	//
	// Now the reaction is set we add out detectors
	//
	
	
	//Add some gamma detectors
	for(int i=1;i<4;i++){
		double angle=pi*i/4;
		TRotation rot;
		rot.RotateY(angle);
		TVector3 vec(120*sin(angle),0,120*cos(angle));
		add_square(experiment,vec,rot,80);
		experiment.set_valid_part(experiment.detN()-1,0,0,0,1);
	}
	
	add_S3(experiment,30);
	//set that detector to only be sensitive to ejectile (target like in our case)
	experiment.set_valid_part(experiment.detN()-1,0,1,0,0);
	
	experiment.print_detectables();
	experiment.print_decay();

	//Draw the geometry in several ways
	can_view->cd(1);
	gPad->Divide(2,2);
	gPad->cd(1);
	experiment.draw_exp();
// 	can_view->cd(1);
// 	gPad->cd(2);
// 	experiment.draw_exp(3);
	can_view->cd(1);
	gPad->cd(2);
	experiment.draw_phi(true);//"true" adds colour

	// Next we spend 20 seconds drawing events so that if something is very wrong it can be corrected
	//because this simple geometry has no obstructions we can run "Obstructions false" mode which is faster
	can_view->cd(2);
// 	experiment.draw_hits_3D(2,2,false,0.5,1,true,20);
// draw_hits_3D(projection,hit_multiplicity,obstructions,display_time,refresh_rate,draw_misses,run_time)
experiment.dist_target_KE_beam.Draw();
	app->Run();
	return 0;
	
	// Draw the target and decay things
	// quite important to check things are set and working right
	can_view->cd(1);
	gPad->cd(3);
	experiment.draw_decay_Z(1,false);
	
	can_view->cd(1);
	gPad->cd(4);
// 	experiment.draw_target_interaction(1,false);
	experiment.draw_target_interaction();
	
	// Next we simulate some events which will be recorded as good and saved to detector
	// only for a ejectile and gamma hit i.e. S3 germanium hit
	experiment.basic_hit_count(1000000,false,0,1,0,1);
	
	// Finally draw detector energy spectra
	// We add some intrinsic resolution to the gamma data
	can_view->cd(2);
	gPad->Clear();
	gPad->Divide(2,2);
	gPad->cd(1);
	experiment.get_det(4).energy.DrawCopy();
	can_view->cd(2);
	gPad->cd(2);
	detector a=experiment.get_det(0);
	add_resolution(&a.energyG)->DrawCopy();
	can_view->cd(2);
	gPad->cd(3);
	detector b=experiment.get_det(1);
	add_resolution(&b.energyG)->DrawCopy();
	can_view->cd(2);
	gPad->cd(4);
	detector c=experiment.get_det(2);
	add_resolution(&c.energyG)->DrawCopy();


	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	