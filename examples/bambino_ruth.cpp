
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
#include <TGraph.h>
#include <TBrowser.h>

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

//Unused historic data from astar
TGraph P_range();
TGraph P_energy();

int main(int argc, char *argv[]){	
// int argcb;
// char **argvb;	
TApplication *app = new TApplication("app", &argc, argv);
TH1::AddDirectory(kFALSE);//avoid name overwrites but could cause other memory histogram issues
//        h->SetDirectory(0);          for the current histogram h


	gStyle->SetOptStat(0);

	double worldsize=80;
	exp_core expr(worldsize);

	add_S3(expr,20,true);
    
	expr.set_beam(66,154);
    expr.set_targ(28,60);
	expr.set_E_beam(3.25*154);
    
    expr.auto_rutherford(10000,true);
    
	expr.set_targ(6,12);
	expr.set_E_beam(3.9*154);
    
    expr.auto_rutherford(10000,true);
    
	expr.set_targ(target(28,60,1));
	expr.set_E_beam(3.25*154);
    
    expr.auto_rutherford(10000,true);
    
	expr.set_targ(target(6,12,1));
	expr.set_E_beam(3.9*154);
    
    expr.auto_rutherford(10000,true);
    


	
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
