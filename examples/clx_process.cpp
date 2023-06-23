
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


bool mask_check(TH1* mask,double& max){
	int X=mask->FindBin(max);
	for(int i=X;i<mask->GetNbinsX();i++){
		if(mask->GetBinContent(i)>0)return true;
	}
	return false;
}



int main(int argc, char *argv[]){	
// int argcb;
// char **argvb;	
TApplication *app = new TApplication("app", &argc, argv);
TH1::AddDirectory(kFALSE);//avoid name overwrites but could cause other memory histogram issues
//        h->SetDirectory(0);          for the current histogram h
	
	int b_z=58,b_a=160;
	target tharget(82,208,2.0);//lead


	TFile infile("inputs/Ex_Ertest.root");
	TFile outfile("outputs/Ex_Ertestout.root","RECREATE");
// 	TFile infile("inputs/Ex_Kr_by_Pd_b78_t110.root");
// 	TFile outfile("outputs/Ex_Kr_by_Pd.root","RECREATE");
	gROOT->cd();
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// END OF USER INPUT ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TCanvas bob;
	exp_core experiment(80.0);

	experiment.set_targ(tharget);

	experiment.set_beam(b_z,b_a,100);
	
	TGraph2D ruth;
	ruth.SetTitle("Ruth");	
	
	TGraph2D ruth_p,ruth_n;
	ruth_p.SetTitle("Rutherford_ped");
	ruth_n.SetTitle("Rutherford_ned");

	for(int e=10 ;e<2000;e+=50){//LOOP ENERGY
		experiment.set_E_beam(e);

		cout<<e<<endl;
		for(double z=10;z<40.0;z+=0.5){//LOOP S3 POSITIONS
			double lower=atan(11.0/z);
			double upper=atan(35.0/z);
			double cxr=experiment.get_rutherford_lab(lower,upper);
			cxr/=1000.0;//mbarn->barn

			ruth_n.SetPoint(ruth_n.GetN(),e+0.001,z+0.001,cxr);

			if(upper>atan(9.0/7.5)&&lower<atan(9.0/7.5)){
				upper=atan(9.0/7.5);
				cxr=experiment.get_rutherford_lab(lower,upper);
				cxr/=1000.0;//mbarn->barn
			}
							
			if(lower<atan(9.0/7.5)){
				ruth_p.SetPoint(ruth_p.GetN(),e+0.001,z+0.001,cxr);
			}	
		}
	}
	
	outfile.cd();
		ruth_p.Draw("tri");
		ruth_p.GetHistogram()->Write(ruth_p.GetTitle());
		ruth_n.Draw("tri");
		ruth_n.GetHistogram()->Write(ruth_n.GetTitle());
	gROOT->cd();	

	
	//set thick target
	experiment.set_target_interaction(2);	
	
	for(int i=0;i<10;i++){//LOOP STATES
		stringstream ss;
		ss<<i;
		string load="Level_"+ss.str();
		if(infile.GetDirectory(load.c_str())){
			cout<<endl<<load;	
			TGraph2D dataa,datab,datac,datad;
			
			dataa.SetTitle((load+"_beam_clx_ped").c_str());
			datab.SetTitle((load+"_beam_clx_ned").c_str());
			datac.SetTitle((load+"_beam_clx_ped_both").c_str());
			datad.SetTitle((load+"_beam_clx_ned_both").c_str());
			
			TGraph2D ratioa,ratiob,ratioc,ratiod;

			ratioa.SetTitle((load+"_clxRuth_ratio_beam_ped").c_str());
			ratiob.SetTitle((load+"_clxRuth_ratio_both_ped").c_str());	
			ratioc.SetTitle((load+"_clxRuth_ratio_beam_ned").c_str());
			ratiod.SetTitle((load+"_clxRuth_ratio_both_ned").c_str());	
						
			for(int e=1;e<10000;e++){//LOOP ENERGY
				stringstream load_g;
				load_g<<load<<"/"<<e<<"_MeV";
				if(infile.Get(load_g.str().c_str())){
					TGraph* clxdist = (TGraph*) infile.Get(load_g.str().c_str());
					clxdist->Draw("al");
					gPad->Update();
					TF1 f1("f",[&](double *x, double *){ return clxdist->Eval(x[0]); },0,180,0);	
					//Because integrals of TGraphs are weird we reference it into a TF1
					double acc=1.e-6;//1.e-12 have to set lower because of unclear error
					double poptot=f1.Integral(0,180,acc);				
						
					experiment.set_E_beam(e);
					experiment.set_primary_dist(clxdist);
					
					cout<<endl<<load_g.str();
					
					double safecm=safe_coulex_angle(b_a,b_z,tharget.A(),tharget.Z(),e); //[Inputs(A1,Z1,A2,Z2,E_lab)]	
					
					for(double z=10;z<40.0;z+=0.5){//LOOP S3 POSITIONS
						
						double lower=atan(11.0/z);
						double upper=atan(35.0/z);
						
						experiment.mask_manual(lower,upper,false,true);
						double pop=poptot*experiment.frac_masked();
						if(mask_check(&(experiment.manual_primary_dist_mask),safecm))pop=0;
						
						experiment.mask_manual(lower,upper,true,true);
						double pop_both=poptot*experiment.frac_masked();
						if(mask_check(&(experiment.manual_primary_dist_mask),safecm))pop_both=0;
						
						double rn=ruth_n.Interpolate(e,z);
						
						datab.SetPoint(datab.GetN(),e,z,pop);
						datad.SetPoint(datad.GetN(),e,z,pop_both);
						if(rn>0){
							ratioc.SetPoint(ratioc.GetN(),e,z,pop*pop/rn);
							ratiod.SetPoint(ratiod.GetN(),e,z,pop_both*pop_both/rn);
						}
						
						
						if(upper>atan(9.0/7.5)&&lower<atan(9.0/7.5)){
							upper=atan(9.0/7.5);
							experiment.mask_manual(lower,upper,false,true);
							pop=poptot*experiment.frac_masked();
							if(mask_check(&(experiment.manual_primary_dist_mask),safecm))pop=0;
							experiment.mask_manual(lower,upper,true,true);
							pop_both=poptot*experiment.frac_masked();	
							if(mask_check(&(experiment.manual_primary_dist_mask),safecm))pop_both=0;	
						}
										
						if(lower<atan(9.0/7.5)){
							rn=ruth_p.Interpolate(e,z);
				
							dataa.SetPoint(dataa.GetN(),e,z,pop);
							datac.SetPoint(datac.GetN(),e,z,pop_both);
							if(rn>0){
								ratioa.SetPoint(ratioa.GetN(),e,z,pop*pop/rn);
								ratiob.SetPoint(ratiob.GetN(),e,z,pop_both*pop_both/rn);
							}
						}
					}
				}
			}//LOOP ENERGY
			
			outfile.cd();
				dataa.Draw("tri");
				dataa.GetHistogram()->Write(dataa.GetTitle());
				datab.Draw("tri");
				datab.GetHistogram()->Write(datab.GetTitle());
				datac.Draw("tri");
				datac.GetHistogram()->Write(dataa.GetTitle());
				datad.Draw("tri");
				datad.GetHistogram()->Write(datab.GetTitle());
				ratioa.Draw("tri");
				ratioa.GetHistogram()->Write(ratioa.GetTitle());
				ratiob.Draw("tri");
				ratiob.GetHistogram()->Write(ratiob.GetTitle());
				ratioc.Draw("tri");
				ratioc.GetHistogram()->Write(ratioc.GetTitle());
				ratiod.Draw("tri");
				ratiod.GetHistogram()->Write(ratiod.GetTitle());

			gROOT->cd();			

		}
	}//LOOP STATES
	
	





	outfile.Close();
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	

	
