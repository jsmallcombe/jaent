
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

	
		
	
	cout<<endl<<"//////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	