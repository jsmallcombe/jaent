
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

#include <james_hist_formatting.h>

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

double E[131]={
0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,17,18,20,22.5,25,27.5,30,32.5,35,37.5,40,45,50,55,60,65,70,80,90,100,110,120,130,140,150,160,170,180,200,225,250,275,300,325,350,375,400,450,500,550,600,650,700,800,900,1000};

double p[131]={0.1342,0.1452,0.1558,0.1662,0.1762,0.1861,0.1957,0.2051,0.2143,0.2323,0.254,0.2751,0.2956,0.3157,0.3355,0.355,0.3744,0.3937,0.432,0.4701,0.5084,0.5468,0.5855,0.6246,0.7041,0.7856,0.8692,0.955,1.04,1.13,1.23,1.32,1.42,1.52,1.62,1.83,2.11,2.4,2.7,3.02,3.35,3.69,4.05,4.41,5.18,5.99,6.85,7.75,8.69,9.67,11.74,13.97,16.33,18.84,21.48,24.28,27.22,30.29,33.5,36.85,40.33,47.69,57.61,68.32,79.8,92.05,105.03,118.76,133.2,148.36,180.77,215.93,253.8,294.32,337.45,383.14,482.01,590.78,709.23,837.16,974.42,1120,1280,1440,1610,1800,1990,2390,2950,3550,4200,4910,5660,6460,7300,8190,10100,12180,14430,16850,19420,22160,28070,34560,41620,49210,57320,65930,75010,84550,94540,104960,115790,138630,169270,202080,236900,273570,311960,351950,393430,436300,525800,619840,717860,819400,924020,1030000,1250000,1480000,1720000};

double d[131]={0.1631,0.1771,0.1907,0.2041,0.2173,0.2303,0.243,0.2555,0.2678,0.2918,0.3206,0.3484,0.3751,0.401,0.426,0.4503,0.474,0.4971,0.542,0.5854,0.6276,0.6689,0.7095,0.7495,0.8285,0.9066,0.9843,1.06,1.14,1.22,1.3,1.38,1.46,1.54,1.62,1.79,2.01,2.24,2.47,2.71,2.95,3.21,3.47,3.74,4.29,4.88,5.49,6.13,6.79,7.48,8.93,10.48,12.11,13.84,15.65,17.54,19.51,21.56,23.68,25.88,28.15,32.91,39.25,46.04,53.29,60.97,69.09,77.63,86.59,95.96,115.9,137.43,160.51,185.12,211.23,238.82,298.28,363.42,434.1,510.21,591.65,678.34,770.18,867.12,969.08,1080,1190,1430,1750,2100,2490,2900,3330,3800,4290,4810,5920,7140,8450,9860,11370,12970,16450,20290,24470,29000,33850,39030,44510,50310,56400,62790,69460,83630,102850,123680,146030,169840,195040,221570,249390,278430,339940,405790,475660,549250,626280,706520,875600,1060000,1240000};

double t[131]={0.1812,0.1972,0.2129,0.2282,0.2433,0.2581,0.2727,0.287,0.3012,0.329,0.3626,0.3953,0.4269,0.4577,0.4875,0.5165,0.5447,0.5721,0.6251,0.6759,0.7246,0.7718,0.8176,0.8623,0.9489,1.03,1.11,1.19,1.27,1.35,1.43,1.51,1.59,1.67,1.74,1.9,2.1,2.31,2.52,2.73,2.95,3.17,3.4,3.63,4.11,4.6,5.12,5.66,6.21,6.79,7.99,9.26,10.6,12.01,13.48,15.02,16.61,18.27,19.99,21.76,23.59,27.41,32.48,37.87,43.57,49.57,55.85,62.44,69.34,76.54,91.81,108.24,125.8,144.47,164.22,185.04,229.8,278.68,331.57,388.4,449.08,513.55,581.75,653.64,729.15,808.25,890.9,1070,1310,1570,1850,2150,2470,2810,3170,3550,4360,5250,6210,7240,8340,9510,12050,14850,17910,21210,24770,28560,32590,36850,41340,46050,50970,61460,75730,91260,107990,125890,144910,165020,186180,208370,255630,306630,361140,418970,479940,543890,679970,826250,981750};

double He4[131]={0.104,0.1141,0.124,0.1338,0.1436,0.1532,0.1627,0.1721,0.1813,0.1994,0.2212,0.2422,0.2625,0.282,0.301,0.3193,0.3372,0.3545,0.3878,0.4194,0.4497,0.4788,0.5067,0.5336,0.585,0.6335,0.6795,0.7236,0.7658,0.8066,0.8462,0.8845,0.9219,0.9584,0.9942,1.06,1.15,1.23,1.31,1.39,1.46,1.54,1.61,1.69,1.83,1.98,2.13,2.27,2.42,2.57,2.88,3.19,3.51,3.84,4.18,4.53,4.89,5.26,5.65,6.04,6.44,7.27,8.37,9.53,10.75,12.04,13.38,14.78,16.24,17.76,20.96,24.38,28.02,31.86,35.9,40.15,49.21,59.05,69.69,81.14,93.36,106.32,120.02,134.44,149.59,165.44,181.99,217.15,264.89,316.78,372.75,432.7,496.59,564.35,635.93,711.28,873,1050,1240,1440,1660,1900,2400,2960,3560,4220,4930,5680,6480,7330,8220,9160,10140,12240,15090,18200,21560,25160,29000,33060,37340,41840,51470,61890,73080,85000,97620,110900,139340,170140,203110};

double Li6[131]={0.0682,0.0748,0.0814,0.0879,0.0944,0.1008,0.1072,0.1136,0.1199,0.1324,0.1479,0.1632,0.1784,0.1934,0.2082,0.2228,0.2371,0.2511,0.2785,0.3048,0.3302,0.3548,0.3787,0.4019,0.4466,0.4892,0.5301,0.5694,0.6073,0.644,0.6796,0.7142,0.7479,0.7808,0.8129,0.8751,0.9496,1.02,1.09,1.16,1.22,1.28,1.34,1.4,1.52,1.63,1.74,1.85,1.95,2.05,2.26,2.45,2.64,2.84,3.03,3.22,3.4,3.59,3.78,3.97,4.16,4.55,5.04,5.54,6.06,6.58,7.12,7.67,8.24,8.81,10.02,11.28,12.59,13.97,15.4,16.9,20.05,23.44,27.06,30.91,34.97,39.26,43.78,48.55,53.56,58.79,64.25,75.81,91.48,108.48,126.77,146.33,167.14,189.19,212.46,236.92,289.37,346.46,408.12,474.25,544.8,619.68,782.09,961.24,1160,1370,1600,1840,2100,2370,2660,2960,3280,3950,4870,5880,6970,8140,9390,10710,12110,13580,16740,20180,23890,27850,32060,36520,46110,56590,67890};

double Li7[131]={0.0691,0.0759,0.0826,0.0893,0.096,0.1027,0.1093,0.1159,0.1225,0.1355,0.1515,0.1673,0.183,0.1986,0.214,0.2292,0.2443,0.2592,0.2882,0.3163,0.3435,0.3699,0.3955,0.4204,0.4683,0.514,0.5578,0.5999,0.6405,0.6798,0.718,0.7551,0.7911,0.8263,0.8607,0.9273,1.01,1.08,1.16,1.23,1.3,1.36,1.43,1.49,1.61,1.73,1.84,1.96,2.06,2.17,2.38,2.58,2.78,2.97,3.17,3.36,3.55,3.74,3.93,4.12,4.3,4.69,5.17,5.65,6.15,6.65,7.17,7.7,8.23,8.78,9.91,11.09,12.32,13.6,14.93,16.31,19.22,22.32,25.63,29.13,32.82,36.71,40.78,45.05,49.51,54.19,59.07,69.42,83.43,98.6,114.91,132.35,150.88,170.51,191.2,212.96,259.55,310.24,364.95,423.6,486.12,552.47,696.3,854.86,1030,1210,1420,1630,1860,2100,2350,2620,2900,3500,4310,5200,6170,7200,8310,9480,10720,12030,14830,17880,21180,24710,28460,32440,41020,50410,60560};

// 	double Si_density=2.3290;//g/cm 3
	double dE_um=150;
	double E_um=1000;
	string beamZ="Ar";
	int beamA=36;
	double stopper_mg=giveme_areal(19.32,25,false)/1000.0;
	double amu=5.2;
	double MeV=120;

	
	//Process command line inputs
	for(int i=1;i<argc;i++){
		string str=argv[i];
		stringstream ss;
		
		//Data file loading
		if(str.find("_")<str.size()){//Beam
			ss<<str.substr(0,str.find("_"));
			ss>>beamA;	
			beamZ=str.substr(str.find("_")+1);
		}else if(str.find("dE")<str.size()){//Set the E and dE thicknesses
			ss<<str.substr(str.find("dE")+2,str.size());
			ss>>dE_um;
		}else if(str.find("E")<str.size()){
			ss<<str.substr(str.find("E")+1,str.size());
			ss>>E_um;
		}else if(str.find("MeV")<str.size()){
			ss<<str.substr(0,str.size()-3);
			ss>>MeV;
		}else if(str.find("amu")<str.size()){
			ss<<str.substr(0,str.size()-3);
			ss>>amu;
		}else{//Stopper
			ss<<str;
			ss>>stopper_mg;
		}
	}

	if(MeV>0)amu=MeV/beamA;

	//
	// Basic exp_core setup
	//

	double worldsize=80;
	exp_core experiment(worldsize);
	TCanvas * can_view = new TCanvas("can_view", "can_view", 1200, 600);
	gStyle->SetOptStat(0);
	can_view->cd();		
	TRandom2 rand;
	rand.SetSeed();
	
	bool random=true;
	bool realtarget=true;
	
	//
	// Set the target
	//
// 	target tharget(1,2,1.0,TVector3(0,0,-1),2,stopper_mg,79,197); //compound 2 polyethelyn
	target tharget(20,40,0.5,TVector3(0,0,-1),0,stopper_mg,79,197);
	
	if(realtarget)experiment.set_targ(tharget);
	else experiment.set_targ(tharget.targ_Z,tharget.targ_A);
	
	experiment.print_target();

	//
	// Set the beam
	//
	
	experiment.set_beam(beamZ,beamA,beamA*amu);
	if(realtarget)experiment.set_target_interaction(2);
	
	//
	// Create output file
	//
    
	stringstream name;
	name<<dE_um<<"um"<<E_um<<"um"<<beamA<<beamZ<<(int)amu*beamA<<"MeV";
	if(stopper_mg>0)name<<stopper_mg<<"mg"<<tharget.backing_A<<nuclear_data_ob::get_symb(tharget.backing_Z);
	name<<".root";
	
	TFile* outfile= new TFile(name.str().c_str(),"RECREATE");
	
	//
	// Create output histograms
	//
		
	outfile->cd();
		TH2F hEtot("Total_Ejectile_Energy","Total_Particle_Energy",1000,0,pi/2,1000,0,50);axislab(&hEtot,"Lab #theta [rad.]","Energy [MeV]");
		TH2F hE("Energy_in_pad","Energy_in_pad",1000,0,pi/2,1000,0,50);axislab(&hE,"Lab #theta [rad.]","Energy [MeV]");
		TH2F hdE("Energy_in_#deltaE","Energy_in_#deltaE",1000,0,pi/2,1000,0,50);axislab(&hdE,"Lab #theta [rad.]","Energy [MeV]");
		TH2F hEpunch("Energy_After_Puncthrough","Energy_After_Puncthrough",1000,0,pi/2,1000,0,50);axislab(&hEpunch,"Lab #theta [rad.]","Energy [MeV]");
		TH2F hRange("Ejectile_Range_in_Silicon","Ejectile_Range_in_Silicon",1000,0,pi/2,1000,0,5000);axislab(&hRange,"Lab #theta [rad.]","Ejectile range in silicon [#mum]");
			
		TH2F hdEE("#deltaE_vs_E","#deltaE_vs_E",1000,0,50,1000,0,30);axislab(&hdEE,"Pad Energy [MeV]","#deltaE Energy [MeV]");
		TH2F hdEdxE("#deltaE/#deltaX_vs_E_sum","#deltaE/#deltaX_vs_E_sum",1000,0,50,1000,0,30);axislab(&hdEdxE,"Energy Sum [MeV]","#deltaE/#deltaX [arb.]");
		TH2F hdEdxEreal("#deltaE/#deltaX_vs_E_sum_Observed","#deltaE/#deltaX_vs_E_sum_Observed",1000,0,50,1000,0,30);axislab(&hdEdxEreal,"Energy Sum [MeV]","#deltaE/#deltaX [arb.]");
		TH3F hdEdxTheta("#deltaE/#deltaX_E_Theta","#deltaE/#deltaX_E_Theta",250,0,25,150,0,15,100,0,pi/2);axislab(&hdEdxTheta,"Energy Sum [MeV]","#deltaE/#deltaX [arb.]","Lab #theta [rad]");
	gROOT->cd();
	

	//
	// Detector Setup
	//

	double ZZ=24.6+7.5;//24.6+7.5 spice new standards
	
	double ta=atan(8./ZZ);
	double tb=atan(38./ZZ);
	
	add_S3(experiment,ZZ,true);
	int last=experiment.detN();
	add_S3(experiment,ZZ+4,false);
	experiment.set_valid_part(last,0,1,0,0);
	for(int i=1;i<last;i++){
		experiment.set_pair(i,last);
	}	
	
	experiment.print_detectables();
	experiment.print_doubles();
	
	can_view->cd();	
	experiment.draw_exp();	
	outfile->cd();
		can_view->Write("Detectors_Geometry");	
	gROOT->cd();
	
	
	//
	//	Setup first reaction
	//
	
	TGraph Si_range(131,E,p);
	TGraph Si_energy(131,p,E);
	experiment.set_ejec(1,1);
	double freeKE=experiment.get_KE_1_tot_CoM();
// 	experiment.set_primary_dist(new TF1("un","1",0,pi));//have to set manual to use mask
// 	experiment.mask_manual(ta,tb,false,true);//mask it to be useful angles
	experiment.print_reaction();	
	can_view->cd();	
	can_view->Clear();	
	experiment.draw_primary_kinematics();
	outfile->cd();
		can_view->Write("KinematicsP");	
	gROOT->cd();

	//
	//	Begin Generating data
	//
	
	int counts=160000;
	for(int i=0;i<counts;i++){
		experiment.basic_make_single(0,0,2);//require dE and E have ejectile//no obstructions
		TLorentzVector ej=experiment.GetLor(1);

		experiment.set_E_star(rand.Uniform(0,freeKE));

		double theta=abs(ej.Theta());
		double thetamid=theta;//Only needed for random addition
		if(random){
			theta+=rand.Gaus(0,0.28);//16 degree angular straggling on exit of target
			thetamid=theta+rand.Gaus(-0.005,0.28);//16 degree angular straggling on exit of dE detector
		}
		
		double E_total=get_KE(&ej);
		double R_total=Si_range.Eval(E_total);
		double E_punch=0;
		double dR=R_total-(dE_um/cos(theta))-(E_um/cos(thetamid));
		if(dR>0)E_punch=Si_energy.Eval(dR);
		double ddR=R_total-((dE_um)/cos(theta));
		double postdE=0;
		if(ddR>0)postdE=Si_energy.Eval(ddR);
		double E_dE=E_total-postdE;
		double E_E=postdE-E_punch;

		if(random){
			//50 keV noise + sqrt(E) res 100KeV @ 5MeV
			if(E_dE>0)E_dE+=rand.Gaus(0,0.05+(sqrt(E_dE)*0.045));
			if(E_E>0)E_E+=rand.Gaus(0,0.05+(sqrt(E_dE)*0.045));
		}		
		if(E_dE<0)E_dE=0;
		if(E_E<0)E_E=0;
		
		
		hEtot.Fill(theta,E_total);
		hEpunch.Fill(theta,E_punch);
		hE.Fill(theta,E_E);
		hdE.Fill(theta,E_dE);
		hRange.Fill(theta,R_total);
		hdEE.Fill(E_E,E_dE);
		
		double thetapixel=atan((floor(ZZ*tan(theta))+0.5)/ZZ);
		hdEdxE.Fill(E_E+E_dE,E_dE*cos(thetapixel));
		if(E_E>0.05&&E_dE>0.05){
			hdEdxEreal.Fill(E_E+E_dE,E_dE*cos(thetapixel));
			hdEdxTheta.Fill(E_E+E_dE,E_dE*cos(thetapixel),theta);
		}
		
		if(i==(int)counts/4){
			experiment.set_ejec(1,2);
			experiment.set_E_star(0);
			freeKE=experiment.get_KE_1_tot_CoM();
			Si_range=TGraph(131,E,d);
			Si_energy=TGraph(131,d,E);
// 			experiment.set_primary_dist(new TF1("un","1",0,pi));//have to set manual to use mask
// 			experiment.mask_manual(ta,tb,false,true);//mask it to be useful angles	
			experiment.print_reaction();		
			can_view->cd();		
			can_view->Clear();	
			experiment.draw_primary_kinematics();
			outfile->cd();
			can_view->Write("KinematicsD");	
			gROOT->cd();
		}
		
		if(i==(int)counts*2/4){
			experiment.set_ejec(1,3);
			experiment.set_E_star(0);
			freeKE=experiment.get_KE_1_tot_CoM();
			Si_range=TGraph(131,E,t);
			Si_energy=TGraph(131,t,E);
// 			experiment.set_primary_dist(new TF1("un","1",0,pi));//have to set manual to use mask
// 			experiment.mask_manual(ta,tb,false,true);//mask it to be useful angles	
			experiment.print_reaction();		
			can_view->cd();		
			can_view->Clear();
			experiment.draw_primary_kinematics();
			outfile->cd();
			can_view->Write("KinematicsT");	
			gROOT->cd();
		}
		
				
		if(i==(int)counts*3/4){
			experiment.set_ejec(2,4);
			experiment.set_E_star(0);
			freeKE=experiment.get_KE_1_tot_CoM();
			Si_range=TGraph(131,E,He4);
			Si_energy=TGraph(131,He4,E);
// 			experiment.set_primary_dist(new TF1("un","1",0,pi));//have to set manual to use mask
// 			experiment.mask_manual(ta,tb,false,true);//mask it to be useful angles		
			experiment.print_reaction();		
			can_view->cd();		
			can_view->Clear();
			experiment.draw_primary_kinematics();
			outfile->cd();
			can_view->Write("KinematicsHe");	
			gROOT->cd();
		}
	}
	

	outfile->cd();
	
	hEtot.Write();
	hEpunch.Write();
	hE.Write();
	hdE.Write();
	hRange.Write();
	hdEE.Write();
	hdEdxE.Write();
	hdEdxEreal.Write();
	hdEdxTheta.Write();

	new TBrowser();
// 	outfile->Close();
	
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
