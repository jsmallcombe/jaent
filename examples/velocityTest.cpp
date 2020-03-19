
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
	
	exp_core experiment(150);
    
	//
	// This target is 2.165 um carbon and downstream 14.917 um of gold backing
	//
	//experiment.set_targ(target(6,12,giveme_areal(2.25,2.165),180.0,0,giveme_areal(19.32,14.917),79,197));

	experiment.set_targ(26,54);
	experiment.set_beam(18,36,115);
    
	//we want an excited beam like so have to force it
	//experiment.set_reco(36,86);    
	experiment.set_ejec(0,2); 
    
	// Now the primary reaction is set, print it to make sure we didnt F up.
	experiment.print_reaction();

// 	TCanvas * can_view = new TCanvas("can_view", "can_view", 1200, 600);
// 	gStyle->SetOptStat(0);
// 	can_view->cd();
// 	can_view->Divide(2,1);
	
	cout<<endl<<"Beta "<<experiment.get_beta_CoM()<<endl;	

	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	
