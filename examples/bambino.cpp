
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
	//
	// Basic exp_core setup
	//

	double worldsize=80;
	exp_core expr(worldsize);
	
	//
	// Set the target
	//
// 	target tharget(1,2,1000,TVector3(0,0,-1),2,stopper_mg*1000,79,197); //compound 2 polyethelyn	
// 	expr.set_targ(target(82,208,1.5));
// 	expr.set_targ(target(28,58,0.98));
// 	expr.set_targ(target(13,27,1.));
	expr.set_targ(target(6,12,0.5));

	//
	// Set the beam
	//
	expr.set_beam(70,158);
// 	expr.set_beam(68,158);
// 	expr.set_beam(54,131);
	
	//
	// Physics
	//
	
	double MeV=safe_coulex_beam(expr.get_BA(),expr.get_BZ(),expr.get_TA(),expr.get_TZ(),TMath::Pi());
// 	double MeV=3.9*expr.get_BA();
// 	double MeV=2.6*expr.get_BA();	
	
	expr.set_E_beam(MeV);
	expr.set_elastic();
// 	expr.set_E_star(3);
	
	expr.set_target_interaction(2);// E^2 target interaction
	
	//
	// Create output file
	//
	
	stringstream name;
	name<<"outputs/"<<expr.channelnames(expr.get_TZ(),expr.get_TA())<<"_"<<expr.channelnames(expr.get_BZ(),expr.get_BA())<<"_"<<(int)MeV<<"MeV.root";
	TFile* outfile= new TFile(name.str().c_str(),"RECREATE");
	

	//
	// Detector Setup
	//

	double ZZa=20;
	double ZZb=20;
	
	double angles[] = {atan(10./ZZa),atan(36./ZZa),atan(10./ZZb),atan(36./ZZb)};
	for(int i=0;i<4;i++)if(angles[i]<0)angles[i]+=pi;
	size_t size = sizeof(angles) / sizeof(angles[0]); 
	sort(angles, angles + size);
	
	add_S3(expr,ZZa,true);
// 	add_S3(expr,ZZb,true);
	
// 	expr.print_detectables();
// 	expr.print_doubles();
	
	
	TCanvas * can_view = new TCanvas("can_view", "can_view", 1200, 600);
	can_view->cd();
		expr.draw_exp();	
	outfile->cd();
		can_view->Write("Detectors_Geometry");	
	gROOT->cd();
	
	//
	//	Setup first reaction
	//
	
	expr.set_primary_dist(new TF1("un","1",0,pi));//have to set manual to use mask
	expr.mask_manual(angles[0],angles[3],true,true);//mask it to be useful angles
	//Not really needed, does save a little on computation especially those small angles which are compressed in CM
	
	expr.print_reaction();	
	
	can_view->cd();	
	can_view->Clear();	
	expr.draw_primary_kinematics();
	outfile->cd();
		can_view->Write("KinematicsP");	
		expr.dist_target_KE_beam.Write("dist_target_KE_beam");
	gROOT->cd();

	//
	//	Begin Generating data
	//
	
// 	new TCanvas;
// 	expr.draw_hits_3D(2,0,true,0.75,1,true,120);
	
	int counts=1000000;
	
	expr.basic_hit_count(counts,true);
	
	//
	// Create output histograms
	//
		
	int N=expr.detN();
	
	outfile->cd();
		TH2F DetE("DetectorsE","DetectorsE",1000,0,MeV,N,0,N);axislab(&DetE,"Energy [MeV]","Detector Number");
		TH2F ETheta("ETheta","ETheta",1000,0,MeV*800,180,0,pi);axislab(&DetE,"Energy [keV]","Theta");
	gROOT->cd();
	
	cout<<endl<<" Reformatting and combining histograms "<<endl;

	
	for(int d=0;d<N;d++){
		detector D=expr.get_det(d);
		
		TH1* DE=add_resolution(&D.energy,3,300);//3% resoultion at 300MeV
		for(int b=1;b<=DE->GetNbinsX();b++){
			double e=DE->GetXaxis()->GetBinCenter(b);
			double n=DE->GetBinContent(b);
			if(n>0){
				DetE.Fill(e,d,n);
			}
		}
		delete DE;
		
		
		TH2* DT=&D.energtTheta;
		
		for(int B=1;B<=DT->GetNbinsY();B++){
			double t=DT->GetYaxis()->GetBinCenter(B);
			
			TH1D* proj=DT->ProjectionX("proj",B,B);
			TH1D* projres = add_resolution(proj,3,300000);
			delete proj;
			
			for(int b=1;b<=projres->GetNbinsX();b++){
				double e=projres->GetXaxis()->GetBinCenter(b);
				double n=projres->GetBinContent(b);
				if(n>0){
					ETheta.Fill(e,t,n);
				}
			}
			delete projres;
		}
	}
	
	outfile->cd();
		DetE.Write();
		ETheta.Write();
	gROOT->cd();
	
	

	outfile->cd();

	new TBrowser();
// 	outfile->Close();
	
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
