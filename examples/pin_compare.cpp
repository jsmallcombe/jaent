
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
	
	gStyle->SetOptStat(0);
	exp_core experiment(60);
	cout<<endl<<"START"<<endl;	
	int count=0;
	
	TFile* output_file= new TFile("pin_compare.root","RECREATE");
	gROOT->cd();

	bool simulate=true;
	
	//
	// Everything above here is basic setup
	//
	

//////////////////////////////////////////////////////
/////////////////  S3   //////////////////////
///////////////////////////////////////////////////////

	spice_auto_setup(experiment,0,true);

	TCanvas * can_view = new TCanvas("S3_view", "S3_view", 1000, 600);
	can_view->Divide(2,1);
	can_view->cd(1);
	experiment.draw_exp();
	can_view->cd(2);
	experiment.draw_exp(3);
	output_file->cd();
		can_view->Write();
	gROOT->cd();
	
	if(simulate){
		TH1D h1 = experiment.theta_cover_auto(true);
		h1.SetName("S3_theta");
		h1.SetTitle("S3_theta");
		h1.SetLineColor(2);
		output_file->cd();
			h1.Write();
		gROOT->cd();
	}	



//////////////////////////////////////////////////////
/////////////////  PIN SETTINGS   //////////////////////
///////////////////////////////////////////////////////


	
	double safety_gap=0.5;

	double A_spacing=10.6+safety_gap;
	double B_spacing=12.4+safety_gap;
	double B_start_offset=12.4+5.3+safety_gap;
	
	double centre_pixel_shift=5.5;

	double Z_nom=28;

	
//////////////////////////////////////////////////////
/////////////////  SFU   //////////////////////
///////////////////////////////////////////////////////

	spice_auto_setup(experiment,0,false);count=0;
	for(int i=-3;i<4;i++){//axis A
		double X=0+A_spacing*i;
		
			for(int j=-3;j<4;j++){//axis B
				
				if(abs(i)==3 && abs(j)==3)continue;
				if(i==0 && j==0)continue;
				
				double Y=B_spacing*j;				
				
				double ys=0,xs=0;
				if(j==0)xs=centre_pixel_shift*abs(i)/i;
				if(i==0)ys=centre_pixel_shift*abs(j)/j;
							
				TVector3 vec(X+xs,Y+ys,Z_nom);
				
				add_PIN(experiment,vec,TRotation(),true);
				count++;
			}			
	}	


	can_view->SetName("SFU_view");
	can_view->SetTitle("SFU_view");
	can_view->cd(1);
	experiment.draw_exp();
	can_view->cd(2);
	experiment.draw_exp(3);
	
	output_file->cd();
		can_view->Write();
	gROOT->cd();
	
	if(simulate){
		TH1D h1 = experiment.theta_cover_auto(true);
		h1.SetName("SFU_theta");
		h1.SetTitle("SFU_theta");
		h1.SetLineColor(1);
		output_file->cd();
			h1.Write();
		gROOT->cd();
	}

	
//////////////////////////////////////////////////////
/////////////////  Quadrants   //////////////////////
///////////////////////////////////////////////////////


output_file->mkdir("Quadrants");
gROOT->cd();
	
string names[8]={"quad_0_","quad_A_","quad_A_1","quad_A_2","quad_A_3","quad_B_","quad_B_1","quad_B_2"};
int A_sh[8]={0,1,1,1,1,0,0,0};
int B_sh[8]={0,0,0,0,0,1,1,1};
int sta[8]={0,-1,1,2,3,-1,1,2};

double Z_step=2;

for(int p=0;p<8;p++){//quad
	//positive for bigger gap. Negative for overlap
	double A_shift=-(1.5+safety_gap)*A_sh[p];//-(1.5+safety_gap); 
	double B_shift=-(2.5+safety_gap)*B_sh[p];//-(2.5+safety_gap); 

	int stack_level=sta[p];
	//Set the level of "shifting" 
	//0 is none
	//-1 is every row/column
	// >0 is various
		
	TRotation rota,rotb;
	rotb.RotateZ(pi);
		
	//Quadrants 
	spice_auto_setup(experiment,0,false);count=0;
	for(int k=0;k<4;k++){//quad
		double Z=Z_nom;
		double X=0;
		
		for(int i=0;i<4;i++){//axis A
			double Y=B_start_offset;
			double dz=0;
			
				for(int j=0;j<3;j++){//axis B
					
					//special rules
					double ys=0;			
					if(i==0)ys=centre_pixel_shift;
					if(i==3 && j==2)continue;
					
					//add detector
					TVector3 vec(-X,Y+ys,Z+dz);
					vec*=rota;
					add_PIN(experiment,vec,rotb,false);
					count++;
					
					//step for next row
					Y+=B_spacing;

					//complex rules to check if next row should have overlap
					//and adjust row and Z accordingly
					if((stack_level==-1)||
					(stack_level==1&&(j>=1||i==3))||	
					(stack_level==2&&j==0)){				
						Y+=+B_shift;
						if(B_shift<-safety_gap){
							dz-=Z_step;
							
						}
					}else{
						dz=0;//If not overlapped got back to nominal
					}
				}
			//step for next column
			X+=A_spacing;

			//complex rules to check if next column should have overlap
			//and adjust column and Z accordingly		
			if((stack_level==-1)||
			(stack_level==1&&i>=1)||
			(stack_level==2&&i<2)||
			(stack_level==3&&i==0)){	
				X+=A_shift;
				if(A_shift<-safety_gap){
					Z-=Z_step;
					
				}
			}else{
				Z=Z_nom;//If not overlapped got back to nominal
			}	
		}
		//Rotate for the next quadrant
		rota.RotateZ(pi/2);
		rotb.RotateZ(pi/2);	
	}
	
	
	////// write /////
	
	can_view->SetName((names[p]+"view").c_str());
	can_view->SetTitle((names[p]+"view").c_str());
	can_view->cd(1);
	experiment.draw_exp();
	can_view->cd(2);
	experiment.draw_exp(3);
	
	output_file->cd("Quadrants");
		can_view->Write();
	gROOT->cd();
	
	
	if(simulate){
		TH1D h1 = experiment.theta_cover_auto(true);
		h1.SetName((names[p]+"theta").c_str());
		h1.SetTitle((names[p]+"theta").c_str());
		h1.SetLineColor(1);
		output_file->cd("Quadrants");
			h1.Write();
		gROOT->cd();
	}
}

	
//////////////////////////////////////////////////////
/////////////////  Quadrants   //////////////////////
///////////////////////////////////////////////////////


output_file->mkdir("Radial");
gROOT->cd();

for(int p=0;p<11;p++){//quad

	Z_nom=20+2*p;

	spice_auto_setup(experiment,0,false);count=0;

	double S3_inner_a=atan(11/(24.6+7.5));
	double pin_inner_r=(Z_nom+7.5)*tan(S3_inner_a);

	double R_base=sqrt((pin_inner_r*pin_inner_r)-((B_spacing/2)*(B_spacing/2)))+(A_spacing/2);
	double r_in=R_base-(A_spacing-safety_gap)/2;

	while((35/(24.6+7.5))>((r_in+0.5)/(Z_nom+7.5))){
		
		double ang=atan(((B_spacing-safety_gap-2.5)/2)/R_base)*2;
		double bang=atan(((B_spacing-safety_gap)/2)/r_in)*2;
		bang-=ang;
		
		int N=(((2*pi)-bang)/ang);
		
		for(int i=0;i<N;i++){//axis A
			if(count>=100)break;
			
			TRotation rotv,rotd;
			rotd.RotateX(-pi/15);
			rotd.RotateZ(pi/2);
			
			rotd.RotateZ(ang*(N-i-1));
			rotv.RotateZ(ang*(N-i-1));
			
			TVector3 vec(0,R_base,Z_nom);
			vec*=rotv;
			add_PIN(experiment,vec,rotd,true);
			count++;

		}
		
		R_base+=(A_spacing+safety_gap+safety_gap);
		r_in=R_base-(A_spacing-safety_gap)/2;
	}
	
	stringstream ss;
	ss<<p;
		
	can_view->SetName(("radial_view"+ss.str()).c_str());
	can_view->SetTitle(("radial_view"+ss.str()).c_str());
	can_view->cd(1);
	experiment.draw_exp();
	can_view->cd(2);
	experiment.draw_exp(3);
	
	output_file->cd("Radial");
		can_view->Write();
	gROOT->cd();
	
	
	if(simulate){
		TH1D h1 = experiment.theta_cover_auto(true);
		h1.SetName(("radial_theta"+ss.str()).c_str());
		h1.SetTitle(("radial_theta"+ss.str()).c_str());
		h1.SetLineColor(1);
		output_file->cd("Radial");
			h1.Write();
		gROOT->cd();
	}
}

		
	output_file->Write();
	output_file->Close();

	cout<<endl<<"////////////////////////////"<<count<<"////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}
