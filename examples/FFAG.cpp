
//
//
// James's not Geant4 : jaent V2.0
// james.smallcombe@outlook.com 01/3/2016
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
	can_view->cd(); 	
	exp_core experiment(250);

	experiment.set_targ(target(nuclear_data_ob::get_Z("Pb"),208,0.166));
// 	experiment.set_targ(6,12);

	cout<<endl<<"Beam";
	experiment.set_beam(1,1,120);
	
	cout<<endl<<"Elastic";
	experiment.set_elastic();
// 	experiment.set_elastic();	
	cout<<endl<<"ElasticDone";

	double angle=pi/4;
	TRotation rot;
	rot.RotateY(angle);
	TVector3 vec(100*sin(angle),0,100*cos(angle));
	add_square(experiment,vec,rot,80);
	experiment.set_valid_part(experiment.detN()-1,0,1,0,0);
	
	experiment.print_reaction();
	experiment.print_detectables();
	experiment.print_target();

// 	experiment.basic_make_single();
// 	experiment.basic_make_single();
// 	experiment.basic_make_single();
// 	experiment.basic_hit_count(100000);
	
	experiment.draw_hits_3D(0,0,true,0.75,1,true,5);

	experiment.auto_rutherford(100000);
// 	experiment.auto_rutherford(int=10000,bool=false);

	cout<<endl<<"//////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	
