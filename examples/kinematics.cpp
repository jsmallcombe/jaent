
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
	
	double worldsize=80;
	exp_core experiment(worldsize);
	
	//
	// Set the target
	//
// 	target tharget(28,58,1.0); //nickel
// 	target tharget(78,196,2.0); //nickel
// 	experiment.set_targ(tharget);

	
	//
	// Set the beam
	//
    

        experiment.set_beam("Dy",154,safe_coulex_beam(154,66,58,28));
        experiment.set_targ("Ni",58);
        experiment.set_elastic();
//         experiment.set_reco("U",241,241.10603330);
//         experiment.set_E_star(0.0);
	
// 	experiment.set_beam("Kr",80,80*4.17);
// 	experiment.set_elastic();
// // 	experiment.set_ejec(1,1);
// 	experiment.set_E_star(1);
	
	experiment.print_reaction();
// 	//
// 	// Input your physics
// 	//
// 	double excited_Z=-1;
// 	double excited_A=-1;
// 	double initial_state_MeV = 1.4;
// 	double transition_keV = 200;
// 	bool excited_beam_like=true;
		// Create the canvas that is used as an inbetween for the whole program
	
	
	TCanvas * can_view = new TCanvas("can_view", "can_view", 1200, 600);
	gStyle->SetOptStat(0);
	can_view->cd();
	
	experiment.draw_primary_kinematics();

	TFile outfile("KGtest.root","RECREATE");
	can_view->Write("kinematics");
	
	
// 	add_S3(experiment,30);
// 	experiment.set_target_interaction(2);
// 	experiment.basic_hit_count(100000);
// 	
// 	for(int i=0;experiment.detN()>i;i++){
// 		stringstream ss;
// 		ss<<"S3element"<<i;
// 		experiment.get_det(i).energy.Write(ss.str().c_str());
// 	}	
	
	outfile.Close();
	
		
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	
