
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

TGraph* clxdis();

double GetDoppler(double E,double beta, TVector3 &ion, TVector3 &gamma) { 
	double tmp = 0;
	double gammma = 1/(sqrt(1-pow(beta,2)));
	tmp = E*gammma *(1 - beta*TMath::Cos(gamma.Angle(ion)));
	return tmp;
}

int main(int argc, char *argv[]){	
// int argcb;
// char **argvb;	
TApplication *app = new TApplication("app", &argc, argv);
TH1::AddDirectory(kFALSE);//avoid name overwrites but could cause other memory histogram issues
//        h->SetDirectory(0);          for the current histogram h



	bool dotemplates=true;//true to generate germanium response template
	bool clxdist=true;
	bool drawcheck=false;//Draw ~20 seconds of check view before full sim
		
	
	//Krypton Experiment
	double MeV=430;//-1 for max safe
	double beamz=36;
	double beama=78;
	target tharget(79,197,1.5);
	vector< double > BeamGammas={};
	vector< double > BGamRelInt={};// should be <= 1000
	vector< double > TargetGammas={547.5,576.5};
	vector< double > TGamRelInt={100,20};// should be <= 1000
	bool detectbeam=true;
	bool spice=true;
	//S3 positions
	double ZZa=32;
	double ZZb=32;
	
	//Erbium Proposal
// 	double MeV=-1;//-1 for max safe
// 	double beamz=68;
// 	double beama=158;
// 	target tharget(28,58,1000);
// 	vector< double > BeamGammas={1454.28,1004.80,1321.2,2775.5};
// 	vector< double > BGamRelInt={10,2,2,1};// should be <= 1000
// 	vector< double > TargetGammas={192.13,335.10,614.26,628.03,820.09,443.13};
// 	vector< double > TGamRelInt={100,20,3,10,5,2};// should be <= 1000
// 	bool detectbeam=false;
// 	//S3 positions
// 	bool spice=false;
// 	double ZZa=32;
// 	double ZZb=20;	
	
	
	gStyle->SetOptStat(0);
	
	//
	// Basic Input
	//

	double worldsize=120;
	exp_core expr(worldsize);
	
	//
	// S3 Detector Setup
	//
	
	double angles[] = {atan(10./ZZa),atan(36./ZZa),atan(10./ZZb),atan(36./ZZb)};
	for(int i=0;i<4;i++)if(angles[i]<0)angles[i]+=pi;
	size_t size = sizeof(angles) / sizeof(angles[0]); 
	sort(angles, angles + size);
	
	add_S3(expr,ZZa,true);
	if(!spice) add_S3(expr,ZZb,true);
	
	int N=expr.detN();
	
	//
	// Use simm to generate theta coverage of the S3 rings
	//
	
	vector< double > elementtheta;
	for(int d=0;d<24;d++){
		elementtheta.push_back(atan((11.5+d)/ZZa));
	}
	for(int d=0;d<24;d++){
		elementtheta.push_back(atan((11.5+d)/ZZb));
	}	
		
	// set the target
	expr.set_targ(tharget);
	// Set the beam
	expr.set_beam(beamz,beama);	
	
	//
	// Physics - Calc the beam energy for safe coulex
	//
	
	if(MeV<0)MeV=safe_coulex_beam(expr.get_BA(),expr.get_BZ(),expr.get_TA(),expr.get_TZ(),TMath::Pi());
	expr.set_E_beam(MeV);
	expr.set_elastic();
	expr.set_target_interaction(2);// E^2 target interaction
	
	
	//Calc kinematics curves (borrowed from exp_core::draw_primary_kinematics()
	double beta = expr.get_beta_CoM();
	TVector3 p(0,0,expr.get_P_1_CoM());
	TGraph angang,anginv,beta1,beta2;
	for(int i=0;i<1001;i++){
		TLorentzVector rl=make_lorentzvec(p,expr.get_BA());
		TLorentzVector el=make_lorentzvec(-p,expr.get_TA());
		el.Boost(0,0,beta);
		rl.Boost(0,0,beta);
		angang.SetPoint(angang.GetN(),rl.Theta(),el.Theta());
		anginv.SetPoint(anginv.GetN(),el.Theta(),rl.Theta());
		beta1.SetPoint(beta1.GetN(),rl.Theta(),get_beta_KE(get_KE(&rl),expr.get_BA()));
		beta2.SetPoint(beta2.GetN(),el.Theta(),get_beta_KE(get_KE(&el),expr.get_TA()));
		p.RotateX(pi/1000.0);
	}
	
	//
	// Calculate mean beta values factors for all the S3 ring angles angles
	//
	
	vector<double> theta_det,theta_bt,theta_tb;
	vector<double> beta_bb,beta_bt,beta_tt,beta_tb;
	for(int d=0;d<N;d++){
		int ring=d;
		if(d>23)ring--;
		
		double dettheta=elementtheta[ring];
		theta_det.push_back(dettheta);
		
		beta_bb.push_back(beta1.Eval(dettheta));//beam-gated beam
		beta_tt.push_back(beta2.Eval(dettheta));//target-gated target
		
		beta_bt.push_back(beta2.Eval(angang.Eval(dettheta)));//beam-gated target
		theta_bt.push_back(angang.Eval(dettheta));
		
		beta_tb.push_back(beta1.Eval(anginv.Eval(dettheta)));//target-gated beam
		theta_tb.push_back(anginv.Eval(dettheta));
	}
	
	//
	// Make a staggered series of Ge Detector response histograms for random generation later
	//
	

	vector< TH1* > gammatemplates;
	if(dotemplates){
		cout<<endl<<"Generating gamma templated";
		for(int i=0;i<40;i++){
			cout<<endl<<"E_gamma "<<((i+1)*50)<<" keV"<<flush;
			stringstream ss;
			ss<<"GTemplate"<<i;
			gammatemplates.push_back(GenGeResponse((i+1)*50));
			gammatemplates[i]->SetName(ss.str().c_str());
		}
	}
	
	//////////////////////////
	////// Build TIGRESS /////
	//////////////////////////
	
	int GeZero=N;
	vector<TVector3> GeCenter;
	vector<double> x={-20,-20,20,20};
	vector<double> y={-20,20,20,-20};
	TVector3 crystal(20.0,20.0,110);
	TVector3 clover(0.0,0.0,110);
	
	for(int j=0;j<4;j++){
		TVector3 place=crystal;
		TVector3 centre=clover;
		TRotation rot;
		rot.RotateY(pi/4);
		rot.RotateZ((j*pi/2)+pi/8);
		place*=rot;
		centre*=rot;
		for(int i=0;i<4;i++){
			GeCenter.push_back(place);
			expr.add_detector(x,y,place,rot,0,0,0,1);
			place.Rotate(pi/2,centre);
			stringstream ss;
			ss<<"Tigress "<<GeCenter.size();
			expr.SetDetName(ss.str());
		}
	}
	
	if(!spice){
		for(int j=0;j<4;j++){
			TVector3 place=crystal;
			TVector3 centre=clover;
			TRotation rot;
			rot.RotateY(3*pi/4);
			rot.RotateZ((j*pi/2)+pi/8);
			place*=rot;
			centre*=rot;
			for(int i=0;i<4;i++){
				GeCenter.push_back(place);
				expr.add_detector(x,y,place,rot,0,0,0,1);
				stringstream ss;
				ss<<"Tigress "<<GeCenter.size();
				expr.SetDetName(ss.str());
				place.Rotate(pi/2,centre);
			}
		}
	}
	
	int dcor=3;
	if(spice)dcor=4;
	
	for(int j=0;j<dcor;j++){
		TVector3 place=crystal;
		TVector3 centre=clover;
		TRotation rot;
		rot.RotateY(pi/2);
		rot.RotateZ((j*pi/4)+pi/8);
		place*=rot;
		centre*=rot;
		for(int i=0;i<4;i++){
			GeCenter.push_back(place);
			expr.add_detector(x,y,place,rot,0,0,0,1);
			stringstream ss;
			ss<<"Tigress "<<GeCenter.size();
			expr.SetDetName(ss.str());
			
			GeCenter.push_back(-place);
			expr.add_detector(x,y,-place,rot,0,0,0,1);
			stringstream SS;
			SS<<"Tigress "<<GeCenter.size();
			expr.SetDetName(SS.str());
			
			place.Rotate(pi/2,centre);
		}
	}

	//////////////////////////
	//////////////////////////
	//////////////////////////
	

	//
	// Create output file and draw the basic things
	//
	
	stringstream name;	name<<"outputs/exbab_"<<expr.channelnames(expr.get_TZ(),expr.get_TA())<<"_"<<expr.channelnames(expr.get_BZ(),expr.get_BA())<<"_"<<(int)MeV<<"MeV.root";
	TFile* outfile= new TFile(name.str().c_str(),"RECREATE");
	
	TCanvas * can_view = new TCanvas("can_view", "can_view", 800, 800);
	
	can_view->cd();	
	can_view->Clear();	
	expr.draw_primary_kinematics();
	outfile->cd();
		can_view->Write("KinematicsP");	
		expr.dist_target_KE_beam.Write("dist_target_KE_beam");
	gROOT->cd();
	
	can_view->cd();	
	can_view->Clear();	
		expr.draw_exp();
	outfile->cd();
		can_view->Write("Geometry");	
	gROOT->cd();
	
	//////////////////////////
	//////////////////////////
	//////////////////////////

	//
	// Begin
	//	

	TF1 stefeffer("stefeff",stefeff,10,3000,4);
	stefeffer.SetParameters(-0.2193,-0.03569,-2434,5.802);
	
	TH1D raw("raw","raw",2000,0,2000);
	TH1D shift("shift","shift",2000,0,2000);
	TH1D shiftO("shiftO","shiftO",2000,0,2000);
	
	
	for(int bt=0;bt<2;bt++){
		bool beamgamma=true;
		bool detectparent=detectbeam;
	
		unsigned int Ny=BeamGammas.size();
		if(Ny>BGamRelInt.size())Ny=BGamRelInt.size();
		
	
		
		if(bt==0){//beam gammas first
			expr.set_ejec(tharget.targ_Z,tharget.targ_A);
		}else{//then targets
			beamgamma=!beamgamma;
			detectparent=!detectparent;
			Ny=TargetGammas.size();
			if(Ny>TGamRelInt.size())Ny=TGamRelInt.size();
			expr.set_ejec(beamz,beama);	
		}
		
		if(clxdist){
			expr.set_primary_dist(clxdis());
			expr.reverse_primary_dist();
		}else  expr.set_primary_dist(new TF1("un","1",0,pi));
		
		expr.mask_manual(angles[0],angles[3],detectparent,!detectparent);//mask it to be useful angles
		//Not really needed, does save a little on computation especially those small angles which are compressed in CM	
			
		//Set S3 ring valids
		for(int d=0;d<N-1;d++){
			expr.set_valid_part(d,0,!detectparent,detectparent,0);
			if(d==23)d++;
		}
		
		int ion=1;
		if(detectparent)ion=2;
		
		
		double intensitysum=0;
		for(unsigned int y=0;y<Ny;y++){
			if(bt==0){
				intensitysum+=BGamRelInt[y];
			}else{
				intensitysum+=TGamRelInt[y];
			}
			
		}
		
		for(unsigned int y=0;y<Ny;y++){
			
			double Ey,Iy;
			if(bt==0){
				Ey=BeamGammas[y];Iy=BGamRelInt[y];
			}else{
				Ey=TargetGammas[y];Iy=TGamRelInt[y];
			}
	
			
			expr.set_E_star(1.5);
			expr.set_gamma(Ey);

			expr.print_reaction();	
			expr.print_decay();
			expr.print_detectables();
			expr.print_target();	

			// TCanvas * viewer = new TCanvas("viewer", "viewer", 800, 800);
			if(drawcheck&&y==0)expr.draw_hits_3D(0,2,true,0.75,1,true,20);
					
			cout<<endl;
			if(beamgamma)cout<<expr.channelnames(expr.get_BZ(),expr.get_BA());
			else cout<<expr.channelnames(expr.get_TZ(),expr.get_TA());
			cout<<" gamma "<<y+1<<"/"<<Ny<<"  "<<Ey<<" keV "<<std::setprecision(2)<< Iy*100/intensitysum<<"% of total"<<endl;
				
				//Screen update
	
			for(int i=0;i<10000*Iy;i++){
				cout << setiosflags(ios::fixed) << i << "/"<<10000*(int)Iy<<"\r"<< flush; 
				expr.basic_make_single(1,0,!detectparent,detectparent,1);//requre gamma and something
					
				vector< vector< int > > det_hits=expr.get_det_hits();
				if(det_hits.size()>3){
					if(det_hits[3].size()>0&&det_hits[ion].size()>0){
						int GeN=det_hits[3][0]-GeZero;
						int ring=det_hits[ion][0];
						if(ring>23)ring--;
						
						if(GeN>=GeCenter.size())continue;
						if(GeN>=elementtheta.size())continue;
						
						TLorentzVector gam=expr.GetLor(3);
						double GammaElab=get_KE(&gam)*1000;
						
						//Generate gamma response function
						double Gelab=GammaElab;
						if(dotemplates){
							int g=round(GammaElab/50);
							if(g<1)g=1;
							if(g>40)g=40;
							g--;
							Gelab=gammatemplates[g]->GetRandom();
							Gelab*=GammaElab/((g+1)*50);
						}
						
						double S3phi=floor(expr.GetLor(ion).Phi()/(pi/16.))*(pi/16.);
						double S3theta=elementtheta[ring];
						
						double ionselfbeta,betanotion,thetanotion;
						if(detectbeam){
							ionselfbeta=beta_bb[ring];
							betanotion=beta_bt[ring];
							thetanotion=theta_bt[ring];
						}else{
							ionselfbeta=beta_tt[ring];
							betanotion=beta_tb[ring];
							thetanotion=theta_tb[ring];
						}
						
						TVector3 ion,notion;
						ion.SetMagThetaPhi(1,S3theta,S3phi);
						notion.SetMagThetaPhi(1,thetanotion,S3phi+pi);
						
		// 				cout<<endl;
		// 				cout<<endl<<ion.Theta()<<" "<<ion.Phi();
		// 				cout<<endl<<expr.GetLor(1).Theta()<<" "<<expr.GetLor(1).Phi();
		// 				cout<<endl<<notion.Theta()<<" "<<notion.Phi();
		// 				cout<<endl<<expr.GetLor(2).Theta()<<" "<<expr.GetLor(2).Phi();

						double GammaEionCM=GetDoppler(Gelab,ionselfbeta, ion, GeCenter[GeN]);
						double GammaEnotionCM=GetDoppler(Gelab,betanotion, notion, GeCenter[GeN]);
						
						raw.Fill(Gelab,stefeffer.Eval(Gelab));
						shift.Fill(GammaEionCM,stefeffer.Eval(Gelab));
						shiftO.Fill(GammaEnotionCM,stefeffer.Eval(Gelab));
					}
				}
			}
		}
	}
	
	
	outfile->cd();
		raw.Write();	
		shift.Write();	
		shiftO.Write();	
	gROOT->cd();

	
// 	



	
	
	


	//
	//	Begin Generating data
	//

	
// 	int counts=1000000;
// 	
// 	expr.basic_hit_count(counts,true);
	
	//
	// Create output histograms
	//
		
	
// 	outfile->cd();
// 		TH2F DetE("DetectorsE","DetectorsE",1000,0,MeV,N,0,N);axislab(&DetE,"Energy [MeV]","Detector Number");
// 	gROOT->cd();
// 	

	
// 	for(int d=0;d<N;d++){
// 		detector D=expr.get_det(d);
// 		
// 		TH1* DE=add_resolution(&D.energy,2,300);//2% resoultion at 300MeV
// 
// 		
// // 		outfile->cd();
// // 		stringstream ss;
// // 		ss<<"Det"<<d;
// // 			DE->Write(ss.str().c_str());
// // 		gROOT->cd();
// // 		
// // 		
// 		for(int b=1;b<=DE->GetNbinsX();b++){
// 			
// 			double e=DE->GetXaxis()->GetBinCenter(b);
// 			double n=DE->GetBinContent(b);
// 			if(n>0){
// 				DetE.Fill(e,d,n);
// 			}
// 		}
// 		
// 		delete DE;
// 	}
// 	
// 	outfile->cd();
// 		DetE.Write();
// 	gROOT->cd();
// 	
// 	
// 
// 	outfile->cd();
// 
	new TBrowser();
// 	outfile->Close();
	
	
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	

TGraph* clxdis(){
	Double_t _fx1[361] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34,34.5,35,35.5,36,36.5,37,37.5,38,38.5,39,39.5,40,40.5,41,41.5,42,42.5,43,43.5,44,44.5,45,45.5,46,46.5,47,47.5,48,48.5,49,49.5,50,50.5,51,51.5,52,52.5,53,53.5,54,54.5,55,55.5,56,56.5,57,57.5,58,58.5,59,59.5,60,60.5,61,61.5,62,62.5,63,63.5,64,64.5,65,65.5,66,66.5,67,67.5,68,68.5,69,69.5,70,70.5,71,71.5,72,72.5,73,73.5,74,74.5,75,75.5,76,76.5,77,77.5,78,78.5,79,79.5,80,80.5,81,81.5,82,82.5,83,83.5,84,84.5,85,85.5,86,86.5,87,87.5,88,88.5,89,89.5,90,90.5,91,91.5,92,92.5,93,93.5,94,94.5,95,95.5,96,96.5,97,97.5,98,98.5,99,99.5,100,100.5,101,101.5,102,102.5,103,103.5,104,104.5,105,105.5,106,106.5,107,107.5,108,108.5,109,109.5,110,110.5,111,111.5,112,112.5,113,113.5,114,114.5,115,115.5,116,116.5,117,117.5,118,118.5,119,119.5,120,120.5,121,121.5,122,122.5,123,123.5,124,124.5,125,125.5,126,126.5,127,127.5,128,128.5,129,129.5,130,130.5,131,131.5,132,132.5,133,133.5,134,134.5,135,135.5,136,136.5,137,137.5,138,138.5,139,139.5,140,140.5,141,141.5,142,142.5,143,143.5,144,144.5,145,145.5,146,146.5,147,147.5,148,148.5,149,149.5,150,150.5,151,151.5,152,152.5,153,153.5,154,154.5,155,155.5,156,156.5,157,157.5,158,158.5,159,159.5,160,160.5,161,161.5,162,162.5,163,163.5,164,164.5,165,165.5,166,166.5,167,167.5,168,168.5,169,169.5,170,170.5,171,171.5,172,172.5,173,173.5,174,174.5,175,175.5,176,176.5,177,177.5,178,178.5,179,179.5,180};
	Double_t _fy1[361] = {0,5.386863e-24,6.660421e-23,2.621043e-18,3.744039e-14,1.112936e-11,4.912762e-10,7.294211e-09,5.503893e-08,2.65042e-07,9.3346e-07,2.621235e-06,6.215562e-06,1.294925e-05,2.437859e-05,4.233937e-05,6.888102e-05,0.0001062211,0.0001566668,0.0002226014,0.0003064277,0.0004105454,0.0005373386,0.0006891743,0.0008683629,0.001077186,0.001317888,0.001592647,0.001903628,0.002252901,0.002642539,0.003074509,0.003550787,0.004073199,0.004643578,0.005263691,0.005935174,0.006659687,0.007438601,0.008273483,0.009165485,0.01011588,0.01112575,0.01219607,0.01332757,0.01452105,0.01577693,0.01709565,0.01847741,0.01992229,0.02143013,0.02300062,0.02463328,0.02632742,0.02808214,0.0298964,0.03176892,0.03369809,0.03568238,0.03771967,0.0398081,0.04194526,0.04412872,0.0463555,0.0486228,0.05092751,0.05326636,0.05563581,0.05803217,0.06045167,0.06289025,0.06534371,0.06780817,0.07027882,0.0727509,0.07522019,0.07768165,0.08013013,0.0825624,0.0849717,0.0873539,0.08970409,0.09201726,0.09428832,0.09651281,0.09868473,0.100801,0.1028553,0.104844,0.1067628,0.1086063,0.1103712,0.1120542,0.1136497,0.1151551,0.1165665,0.1178813,0.1190962,0.120209,0.1212164,0.1221169,0.1229091,0.1235892,0.1241583,0.1246144,0.1249564,0.1251841,0.1252973,0.1252969,0.1251822,0.124955,0.1246148,0.1241647,0.123606,0.1229394,0.1221682,0.1212943,0.1203204,0.1192495,0.1180844,0.1168293,0.1154864,0.1140601,0.1125547,0.1109729,0.1093203,0.1075996,0.1058168,0.1039746,0.1020793,0.1001339,0.09814429,0.0961139,0.09404785,0.0919507,0.08982687,0.08768086,0.08551719,0.08334018,0.08115394,0.07896287,0.07677078,0.07458171,0.07239944,0.07022756,0.06806966,0.06592887,0.06380839,0.06171135,0.05964037,0.05759825,0.05558719,0.05360967,0.05166766,0.049763,0.04789741,0.04607248,0.04428939,0.04254946,0.04085351,0.03920241,0.0375969,0.03603729,0.03452415,0.03305742,0.03163731,0.03026357,0.02893626,0.02765475,0.02641896,0.02522808,0.02408164,0.02297892,0.02191917,0.0209015,0.01992501,0.01898889,0.01809196,0.01723338,0.01641204,0.01562693,0.01487683,0.01416079,0.01347766,0.01282646,0.01220592,0.01161513,0.01105301,0.01051847,0.01001061,0.009528415,0.00907093,0.008637253,0.008226504,0.007837829,0.007470431,0.007123519,0.006796363,0.006488256,0.006198544,0.005926567,0.005671736,0.005433478,0.005211272,0.005004581,0.004812952,0.004635914,0.004473051,0.004323928,0.004188184,0.004065434,0.003955319,0.003857502,0.003771661,0.003697462,0.003634598,0.003582762,0.003541643,0.003510948,0.003490372,0.003479618,0.003478377,0.003486347,0.003503217,0.003528662,0.003562366,0.003603995,0.003653221,0.003709698,0.003773078,0.003842985,0.003919073,0.004000968,0.004088276,0.004180618,0.004277603,0.004378837,0.004483907,0.004592424,0.004703972,0.004818133,0.004934518,0.005052698,0.005172268,0.005292828,0.005413975,0.005535311,0.005656454,0.005777003,0.005896597,0.006014859,0.006131443,0.006246003,0.006358192,0.006467715,0.006574239,0.006677494,0.006777191,0.006873073,0.006964903,0.007052439,0.007135484,0.007213838,0.007287327,0.0073558,0.007419113,0.007477136,0.00752978,0.007576944,0.007618576,0.007654603,0.007685014,0.00770978,0.00772889,0.007742333,0.007750191,0.007752475,0.007749272,0.007740602,0.007726556,0.007707232,0.007682727,0.007653208,0.007618707,0.007579394,0.007535437,0.007486938,0.007434,0.007376921,0.007315737,0.007250613,0.007181776,0.007109386,0.007033539,0.006954455,0.006872346,0.006787361,0.006699592,0.006609297,0.006516563,0.006421638,0.006324601,0.006225693,0.006124985,0.006022654,0.00591882,0.005813649,0.005707301,0.005599836,0.005491451,0.005382199,0.005272221,0.005161609,0.005050493,0.004938912,0.00482699,0.004714843,0.004602499,0.004490026,0.004377489,0.004265024,0.004152574,0.004040265,0.003928113,0.003816195,0.003704488,0.003593065,0.00348193,0.003371148,0.003260689,0.003150581,0.003040852,0.002931523,0.002822561,0.002713986,0.002605838,0.002498059,0.002390661,0.002283656,0.002177029,0.002070768,0.001964875,0.001859317,0.001754096,0.001649195,0.001544597,0.001440284,0.001336239,0.001232457,0.001128905,0.001025576,0.0009224408,0.0008194895,0.0007167011,0.0006140493,0.0005115209,0.
0004090951,0.000306751,0.0002044658,0.0001022234,2.555502e-05};
	TGraph *graph = new TGraph(361,_fx1,_fy1);
	
	return graph;
}