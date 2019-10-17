
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
//        h->SetDirectory(0);          for the current histogram hplace

	double worldsize=120;
	exp_core expr(worldsize);
	
	//Krypton Experiment
// 	double MeV=350;//-1 for max safe
// 	double beamz=28;
// 	double beama=84;
// 	target tharget(20,48,1.0);
    
    double MeV=250;//-1 for max safe
	double beamz=18;
	double beama=36;
	target tharget(20,40,0.1);

	// set the target
	expr.set_targ(tharget);
	expr.set_beam(beamz,beama);	
	expr.set_E_beam(MeV);
	expr.set_elastic();
	expr.set_target_interaction(2);// E^2 target interaction
//     expr.set_ejec(0,4);
    expr.set_ejec(2,4);
	
	//Calc kinematics curves (borrowed from exp_core::draw_primary_kinematics()
	double beta = expr.get_beta_CoM();
	
    vector< TH1* > gammatemplates;
    cout<<endl<<"Generating gamma templated";
    for(int i=0;i<40;i++){
        cout<<endl<<"E_gamma "<<((i+1)*50)<<" keV"<<flush;
        stringstream ss;
        ss<<"GTemplate"<<i;
        gammatemplates.push_back(GenGeResponse((i+1)*50));
        gammatemplates[i]->SetName(ss.str().c_str());
    }
	
	//////////////////////////
	////// Add Detectors /////
	/////////////////////////
    
    vector< TVector3 > dpos;
    TVector3 zero(0,0,1);
    
    TVector3 place(0,0,60);
    TRotation rotP,rot;
    rot.RotateY(-pi/4);
    rotP.RotateY(-pi/4);
    place*=rotP;
    add_annulus(expr,3,15,place,rot);
    dpos.push_back(place);
    rot.RotateY(-pi/4);
    place*=rotP;
    add_annulus(expr,3,15,place,rot);
    dpos.push_back(place);
    rot.RotateY(-pi/4);
    place*=rotP;
    add_annulus(expr,3,15,place,rot);
    dpos.push_back(place);
        
    expr.set_valid_part(0,0,0,0,1);
    expr.set_valid_part(1,0,0,0,1);
    expr.set_valid_part(2,0,0,0,1);

        
	//
	// Create output file and draw the basic things
	//
	TFile* outfile= new TFile("outputs/StAndrews.root","RECREATE");
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

	TF1 stefeffer("stefeff",stefeff,10,3000,4);
	stefeffer.SetParameters(-6.69572e+00,1.04534e+00,-1.23989e+04,1.52477e+10);

	

// 		unsigned int Ny=

// 		if(clxdist){
// 			expr.set_primary_dist(clxdis());
// 			expr.reverse_primary_dist();
// 		}else  expr.set_primary_dist(new TF1("un","1",0,pi));
		

        expr.set_E_star(1.5);

        expr.print_reaction();	
        expr.print_detectables();
        expr.print_target();	

//         vector<double> eng={645.8,784.6,440.3,1224.0,666.7,237.9,765.0,450.4,537.6,68.7};
        vector<double> eng={944.51,1094.4,438.9,964.39,620.7,1438.1,782.6};
        
        for(unsigned int k=0;k<eng.size();k++){
            stringstream ssA;
            ssA<<"A"<<k;
//             TH1D  A;
            stringstream ssB;
            ssB<<"B"<<k;
//             TH1D B
            stringstream ssC;
            ssC<<"C"<<k;
            TH1D HT[3]={TH1D(ssA.str().c_str(),ssA.str().c_str(),2000,0,2000),TH1D(ssB.str().c_str(),ssB.str().c_str(),2000,0,2000),TH1D(ssC.str().c_str(),ssC.str().c_str(),2000,0,2000)};
//             TH1D C;
        
            
            
            expr.set_gamma(eng[k]);
            expr.print_decay();
        
            int Ngam=300000;
			for(int i=0;i<Ngam;i++){
                if(i%100==0)cout<<std::setprecision(2)<< (double)i*100/Ngam<<"%"<<"\r"<<flush;
                
				expr.basic_make_single(1,0,0,0,1);//requre gamma and something
					
				vector< vector< int > > det_hits=expr.get_det_hits();
				if(det_hits.size()>3){
   
                    int id=det_hits[3][0];
					TLorentzVector gam=expr.GetLor(3);
                    double GammaElab=get_KE(&gam)*1000;
                    double Gelab;
                    
                    int g=round(GammaElab/50);
                    if(g<1)g=1;
                    if(g>40)g=40;
                    g--;
                    Gelab=gammatemplates[g]->GetRandom();
                    Gelab*=GammaElab/((g+1)*50);         
                    
                    HT[id].Fill(Gelab);
// 						raw.Fill(GammaElab);
//                         shift.Fill(GammaE);
                    // double GammaE=GetDoppler(GammaElab,beta,zero,dpos[id]);
                    // double GammaEnotionCM=GetDoppler(Gelab,betanotion, notion, GeCenter[GeN]);
                    // raw.Fill(Gelab,stefeffer.Eval(Gelab));
                    // shift.Fill(GammaEionCM,stefeffer.Eval(Gelab));
                    // shiftO.Fill(GammaEnotionCM,stefeffer.Eval(Gelab));
                }
            }
            outfile->cd();
            HT[0].Write();
            HT[1].Write();
            HT[2].Write();
            gROOT->cd();
        }
// 	raw.DrawCopy();
//     shift.SetLineColor(2);
// 	shift.DrawCopy("same");
// 		raw.Write();	
// 		shift.Write();	
// 		shiftO.Write();	
// 	gROOT->cd();

	
// 	


	


	//
	//	Begin Generating data
	//

	
// 	int counts=1000000;
// 	
// 	expr.basic_hit_count(counts,true);

// 	outfile->cd();
// 		TH2F DetE("DetectorsE","DetectorsE",1000,0,MeV,N,0,N);axislab(&DetE,"Energy [MeV]","Detector Number");
// 	gROOT->cd();
    
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
// 	new TBrowser();
// 	outfile->Close();
	
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
// 	app->Run();
	return 0;
}	
