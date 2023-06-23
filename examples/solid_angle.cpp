
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
	
	TCanvas * can_view = new TCanvas("can_view", "can_view", 800, 600);
	gStyle->SetOptStat(0);
	can_view->cd();
	can_view->Divide(3,2);
	exp_core experiment(100);
	
	cout<<endl<<"START"<<endl;	
	
	//
	// Everything above here is basic setup
	//
	
	//As we are only concerned with geometry we create a default gun
	experiment.set_gun();
	//Uniform is default so this is not needed, just speeds example.
	double solid_frac=experiment.set_uniform(0,0.25*pi);
	experiment.reverse_primary_dist();//as we want recoils (see detectors)
	
	
	//
	// This creates a SPICE pin diode detector directly in front of the source 20mm
	// Detector has 2 elements. 0=pcb 1=detector
	//	
	add_PIN(experiment,TVector3(0,0,20),TRotation(),true);
	
	//Make sure the detector is only sensitive to one particle (otherwise we double count)
	experiment.set_valid_part(1,1,0,0,0);
	 
	//Draw the geometry in several ways
	can_view->cd(1);
	experiment.draw_exp();
	can_view->cd(2);
	experiment.draw_exp(3,false);
	can_view->cd(3);
	experiment.draw_phi(true);//"true" adds colour
	
	//
	// Calculate detector coverage
	//
	
	int reps=500000;
	
	// This line calculates the number (and fraction) of events that result in each type of particle being detected
	// in ANY valid detector, and prints it to the terminal
	// This function also adds hit data to detector objects themselves
	experiment.basic_hit_count(reps,false,1);
	// We run in non-obstruction mode for speed
	
	// This line counts the number of events for which detector #specified is hit my any of its valid particles
	double frac_hit=experiment.basic_det_hit_frac(reps,1,false);
	// We run in non-obstruction mode for speed
	
	// Now using the previous run we calculate the solid angle
	double FF=frac_hit*solid_frac;
	double ff=sqrt((1/reps)+(1/(reps*frac_hit)))*FF;
	
	cout<<endl<<endl<<"Solid angle of selected detector is "<<FF*100<<"("<<ff<<") %  "<<frac_hit*solid_frac*4*pi<<" [strad]"<<endl;
		
	// Next we fetch the detector and draw it's hits stored during the the basic_hit_count experiment commands
	can_view->cd(4);
	experiment.get_det(1).hit_pat.DrawCopy("colz");
	can_view->cd(5);
	experiment.get_det(1).energy.DrawCopy();
	
	// We can also get the number of hits from the detector directly, like so
	
	double NN=experiment.get_det(1).energy.Integral();
	
	cout<<endl<<endl<<"Solid angle of selected detector is from histogram "<<NN*solid_frac/reps*100<<" %  "<<NN*solid_frac/reps*4*pi<<" [strad]"<<endl;
	
	// Finally we make the PCB sensitive to particles in order to see its hit pattern
	// But with obstruction ON to demonstrate the effect
	experiment.set_valid_part(0,1,1,0,0);
	experiment.basic_hit_count(reps,true,1);
	can_view->cd(6);
	experiment.get_det(0).hit_pat.DrawCopy("colz");
	

	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	