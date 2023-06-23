
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
	
	
	
// 	TCanvas * can_view = new TCanvas("can_view", "can_view", 800, 600);
// 	can_view->cd(); 	
	exp_core experiment(250);

	
	experiment.set_targ(36,88);		
	experiment.set_beam(6,12,500);
	experiment.set_elastic();
// 	
// 	experiment.print_reaction();	
// 	cout<<endl<<endl<<endl;
// 
// 	
// 	experiment.set_E_star(2.15554);
// 	experiment.set_gamma(3011.45);
// 	experiment.print_decay();	
// 	cout<<endl<<endl<<endl;
// 	
// 	experiment.set_E_star(2.15455);
// 	experiment.set_gamma(215455);	
// 	experiment.print_decay();
// 	cout<<endl<<endl<<endl;	
// 	
// 	experiment.set_E_star(2.15455);
// 	experiment.set_gamma(1225.34);
// 	experiment.print_decay();	
// 	cout<<endl<<endl<<endl;
// 
// 	experiment.set_E_star(2.155322);
// 	experiment.set_ICE(3011.2211);
// 	experiment.print_decay();	
// 	cout<<endl<<endl<<endl;
// 	
// 	experiment.set_E_star(2.1552228);
// 	experiment.set_ICE(2155.232);
// 	experiment.print_decay();	
// 	cout<<endl<<endl<<endl;	
// 	
// 	experiment.set_E_star(2.155198);
// 	experiment.set_ICE(1225.333);
// 	experiment.print_decay();	
// 	cout<<endl<<endl<<endl;
// 	

 	experiment.set_E_star(1.20);	
 	experiment.set_alpha(0.55,1.44);
	experiment.print_decay();	
	
 	experiment.set_E_star(1.20);	
 	experiment.set_b_ray(true,1.5);
	experiment.print_decay();	
	experiment.beta_dist_hist.DrawCopy();
 	experiment.set_E_star(1.20);	
 	experiment.set_b_ray(false,1.5);
	experiment.beta_dist_hist.SetLineColor(2);	
	experiment.beta_dist_hist.DrawCopy("same");
	experiment.print_decay();	
	
	// 	void set_alpha(double=0.0,double=0.0);
	
	
		
	
	cout<<endl<<"//////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	