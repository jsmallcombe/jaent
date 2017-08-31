
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
int main(int argc, char *argv[]){	
// int argcb;
// char **argvb;	
TApplication *app = new TApplication("app", &argc, argv);
TH1::AddDirectory(kFALSE);//avoid name overwrites but could cause other memory histogram issues
//        h->SetDirectory(0);          for the current histogram h
	
	double worldsize=80;
	exp_core experiment(worldsize);
	
	//
	// Set the option you wish to run
	//
	int control =0;
	// 0 = EVERYTHING
	// 1 = geometry
	// 2 = Rutherford
	// 3 = Decay position
	// 4 = Distribution hits
	// 5 = Draw events for real time viewing
	
	bool save_text=true;//saves terminal output but no terminal output until completion
	
	//
	// Set the target
	//

// 	//1.5 mg/cm2 Pt196 + 2 um Pb 207
// 	double density_targ=21.4*(196/195.084);//adjusted for enrichment g/cm3	
// 	double density_back=11.34;//g/cm3
// 	target tharget(78,196,1.5,180.0,0,giveme_areal(11.34,2),82,207);

// //1.6 mg/cm2 Pd110
 	double density_back=0;//g/cm3
	double density_targ=11.9*(110/106.42);//adjusted for enrichment g/cm3
	target tharget(46,110,1600);
	
//1.0 um Au Natural (197)
//  	double density_back=0;//g/cm3
// 	double density_targ=19.32;//adjusted for enrichment g/cm3
// 	target tharget(79,197,giveme_areal(density_targ,1.0));
	
	experiment.set_targ(tharget);

	//
	// Set the beam
	//
	double ebeam=296.4; //MeV
// 	double ebeam=429; //MeV
	experiment.set_beam(36,78,ebeam);	
		
	//
	// Input your physics
	//
	double excited_Z=-1;
	double excited_A=-1;
	double initial_state_MeV = 1.0;
	double transition_keV = 1000;
	bool excited_beam_like=true;
	double halflife=0.0012;//in ns
	string dist="";// leave blank for a uniform distribution
// 	string dist="examples/Kr.txt";// leave blank for a uniform distribution

	//
	// Setup Detectors
	//	
	
	int config=0; // 0 normal, 1 stopper, 2 old
	int S3=2; // 0 off, 1 S3 simple, 2 S3 rings , 3 pin
	double overr=18.5+1+0.85; // mm >0 non-default S3 / pin position
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// END OF USER INPUT ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double t_d=0.2,t_u=0.35*pi;

	///////// First we sort output files ////////
	
	// Create the names 
	stringstream ss;
	ss<<"outputs/SPICE_"<<experiment.channelnames(experiment.get_BZ(),experiment.get_BA());
	ss<<"_"<<experiment.channelnames(experiment.get_TZ(),experiment.get_TA());
	if(tharget.backing_thickness>0)
		ss<<"_"<<experiment.channelnames(tharget.backing_Z,tharget.backing_A)<<"_backed";
	ss<<"_"<<ebeam<<"MeV";
	string file_title=ss.str();
	
	// Check if the root file exists already 
	bool exist=false;
	TFile* outfile= new TFile((file_title+".root").c_str(),"READ");
	if(outfile->IsOpen()){exist=true;outfile->Close();delete outfile;}
	
	// Open the file root file to right, whether or not it exists
	outfile= new TFile((file_title+".root").c_str(),"UPDATE");//Open an existing file for writing. If no file exists, it is created.
	gROOT->cd();// 

// This code is for capturing the terminal output
std::streambuf *psbuf, *backup;
std::ofstream filestr;	
if(save_text){
	// Append if the root file is being
	if(exist){filestr.open(file_title+".txt",std::ofstream::app);filestr<<endl<<"/////// APPEND //////"<<endl;} //appended
	else filestr.open(file_title+".txt");
	backup = std::cout.rdbuf();     // back up cout's streambuf
	psbuf = filestr.rdbuf();        // get file's streambuf
	std::cout.rdbuf(psbuf);         // assign streambuf to cout
}

	// Create the canvas that is used as an inbetween for the whole program
	TCanvas * can_view = new TCanvas("can_view", "can_view", 1200, 600);
	gStyle->SetOptStat(0);
	can_view->cd();
	
	/////////////////////////////////////////////////////
	////////////////// Draw Geometry ////////////////////
	/////////////////////////////////////////////////////
	
	// This loop outputs graphics of the chosen detector geography
	if(control==0||control==1){spice_auto_setup(experiment,config,S3,overr);

		if(!outfile->GetDirectory("Geometry"))outfile->mkdir("Geometry");
		gROOT->cd();
		
		can_view->cd();can_view->Clear();
		experiment.draw_exp();
		outfile->cd("Geometry");
		can_view->SetName("3Dview");can_view->Write("3Dview",TObject::kOverwrite);
		gROOT->cd();
		
		can_view->cd();can_view->Clear();
		experiment.draw_exp(2);
		outfile->cd("Geometry");
		can_view->SetName("side_view");can_view->Write("side_view",TObject::kOverwrite);
		gROOT->cd();
		
		can_view->cd();can_view->Clear();
		experiment.draw_exp(3);
		outfile->cd("Geometry");
		can_view->SetName("beam_view");can_view->Write("beam_view",TObject::kOverwrite);
		gROOT->cd();
		
		can_view->cd();can_view->Clear();
		experiment.draw_phi();
		outfile->cd("Geometry");
		can_view->SetName("theta_view");can_view->Write("theta_view",TObject::kOverwrite);
		gROOT->cd();
	}
	
	/////////////////////////////////////////////////////	
	//////// Calculate Rutherford and Kinematics ////////
	/////////////////////////////////////////////////////

	// This loop calculates the basic beam target kinematics and Rutherford scattering rates
	if(control==0||control==2){spice_auto_setup(experiment,config,S3,overr);
		//experiment.set_rutherford();//sets reaction
		experiment.set_elastic();
		
		can_view->Clear();
		can_view->Divide(2);
		can_view->cd(1);
		experiment.set_ruthmask(t_d,t_u,true,false);// Cut down on wasted simulations by only focusing on the region of interest
		experiment.draw_hit_pattern_2D(1,1000000,1,1,0);
		((TH2D*)gPad->GetPrimitive("Hit_Map_copy"))->SetTitle("Target");
		can_view->cd(2);
		experiment.set_ruthmask(t_d,t_u,false,true);
		experiment.draw_hit_pattern_2D(1,1000000,1,0,1);
		((TH2D*)gPad->GetPrimitive("Hit_Map_copy"))->SetTitle("Beam");
		outfile->cd();
		can_view->SetName("Hit_map_rutherford");can_view->Write("Hit_map_rutherford",TObject::kOverwrite);
		gROOT->cd();
		
		can_view->Clear();can_view->cd();
		experiment.draw_primary_kinematics();
		outfile->cd();
		can_view->SetName("Kinematics");can_view->Write("Kinematics",TObject::kOverwrite);
		gROOT->cd();
		
		experiment.auto_rutherford_OTT(100000,true);
		if(S3==2){spice_auto_setup(experiment,config,1,overr);experiment.auto_rutherford_OTT(100000,true);}
	}
	

	
	/////////////////////////////////////////////////////	
	////// Experimental Interaction and Detection ///////
	/////////////////////////////////////////////////////

	// This loop calculates data for the reaction of interest
	// 
	if(control==0||control==3||control==4||control==5){
		if(excited_Z<0||excited_A<0){
			if(excited_beam_like){
				experiment.set_reco(experiment.get_BZ(),experiment.get_BA());
			}else{
				experiment.set_elastic();
			}			
		}else{
			experiment.set_reco(excited_Z,excited_A);
		}
		
		// Now we set the target interaction and decay information		
		experiment.set_target_interaction(2);//E^2 is a reasonably approximation for most of our stuff
		
		experiment.set_E_star(initial_state_MeV);
		
		experiment.set_ICE(transition_keV);
		// experiment.density_targ=density_targ;// experiment.density_back=density_back;
		if(halflife>0) experiment.set_halflife_ns(halflife,density_targ,density_back); // Can set the densities directly as well

		//set some distribution
		if(dist.size()>0){
			experiment.set_primary_dist(dist);
		}else{
			//experiment.set_uniform(t_d,t_u);
			TFormula fot("fot","sin(x)");
			experiment.set_primary_dist(&fot);//using this one so we can use the same masking for both
		}
		
		if(excited_beam_like)experiment.reverse_primary_dist();//because distribution defined for beam like by default but RECOIL is excited
		
		experiment.set_implant_escape(-1); //SET THIS so things dont escape when they go wack into things
		
		int z=2;// Detector number of the first heavy ion detector
		if(S3==2)z++;// If S3, skip the PCB. PCBs left out for the diodes in spice_auto_setup
		
		//////////////////////////////
		////// FINISHED SET UP ///////
		//////////////////////////////
		
		experiment.print_reaction();
		experiment.print_target();
		experiment.print_decay();
// 		experiment.print_detectables();
		
		outfile->cd();
			experiment.dist_target_KE_beam.Write("KE_though_target",TObject::kOverwrite);
		gROOT->cd();	
		
		if(control==0||control==3){spice_auto_setup(experiment,config,S3,overr);
			experiment.basic_decay_check(100000);

			// Draw the target and decay things
			// quite important to check things are set and working right
			can_view->Clear();can_view->Divide(2);
			can_view->cd(1);
			experiment.draw_decay_Z(0,true,10);
			can_view->cd(2);
			experiment.draw_target_interaction(0,true);
			outfile->cd();
			can_view->SetName("Ungated_Decay_Position");can_view->Write("Ungated_Decay_Position",TObject::kOverwrite);
			gROOT->cd();
			
			for(int i=z;i<experiment.detN();i++){ //set only the daughter recoil as valid in the S3 or pins
				experiment.set_valid_part(i,0,0,1,0);
			}
			string n=experiment.channelnames(experiment.get_D1Z(),experiment.get_D1A())+"_Gated_Decay_Position";
			experiment.mask_manual(t_d,t_u,true,false);//saves on wasted computation

			can_view->Clear();can_view->Divide(2);
			can_view->cd(1);
			experiment.draw_decay_Z(1,true,10);
			can_view->cd(2);
			experiment.draw_target_interaction(1,true);
			outfile->cd();
			can_view->SetName(n.c_str());can_view->Write(n.c_str(),TObject::kOverwrite);
			gROOT->cd();		
			
			
			for(int i=z;i<experiment.detN();i++){ //set only the ejectile as valid in the S3 or pins
				experiment.set_valid_part(i,0,1,0,0);
			}		
			string N=experiment.channelnames(experiment.get_EZ(),experiment.get_EA())+"_Gated_Decay_Position";
			experiment.mask_manual(t_d,t_u,false,true);//saves on wasted computation

			
			can_view->Clear();can_view->Divide(2);
			can_view->cd(1);
			experiment.draw_decay_Z(1,true,10);
			can_view->cd(2);
			experiment.draw_target_interaction(1,true);
			outfile->cd();
			can_view->SetName(N.c_str());can_view->Write(N.c_str(),TObject::kOverwrite);
			gROOT->cd();
		}
		
		// Distribution hit loops
		if(control==0||control==4){spice_auto_setup(experiment,config,S3,overr);
			can_view->Clear();
			can_view->cd();
			
			experiment.mask_manual(t_d,t_u,true,true);
			double ra=experiment.frac_masked();//Ask what fraction we have masked too, as unmask integral = 1	
			// Print a basic output of fractions of hits
			cout<<endl<<endl<<"Multiply by fraction "<<ra;
			experiment.basic_hit_count(100000,true,1,1,0,0);
						
			string r=experiment.channelnames(experiment.get_RZ(),experiment.get_RA());
			string e=experiment.channelnames(experiment.get_EZ(),experiment.get_EA());
			experiment.manual_primary_dist_hist.SetLineWidth(2);			
			
			TH1D full;
			// Get theta CM for different hit conditions					
			// Fetch the manual CM distribution histogram and the speed masked one
			// We want to show the mask so we know we havent cut it too small			
			can_view->Clear();
			can_view->Divide(2);
			can_view->cd(1);
				experiment.mask_manual(t_d,t_u,true,false);//saves on wasted computation
				full=experiment.manual_primary_dist_store;
				ra=experiment.manual_primary_dist_hist.Integral();//Ask what fraction we have masked too, as unmask integral = 1			
				full.SetTitle(("CM_theta_dist_"+r+"_gated").c_str());
				
				//Get theta CM for different hit condition	
				TH1D recoh=experiment.theta_hist(500000,true,1,0,0,0);
				recoh.SetLineColor(2);		
				recoh.Scale(ra*(full.GetBinWidth(1)/recoh.GetBinWidth(1)));//Some basic scaling to match distribution histogram	

				full.DrawCopy();
				experiment.manual_primary_dist_hist.DrawCopy("same");
				recoh.DrawCopy("same");
			can_view->cd(2);
				experiment.mask_manual(t_d,t_u,false,true);//saves on wasted computation
				full=experiment.manual_primary_dist_store;
				ra=experiment.manual_primary_dist_hist.Integral();
				full.SetTitle(("CM_theta_dist_"+e+"_gated").c_str());
				
				TH1D ejech=experiment.theta_hist(500000,true,0,1,0,0);
				ejech.SetLineColor(2);				
				ejech.Scale(ra*(full.GetBinWidth(1)/ejech.GetBinWidth(1)));
				
				full.DrawCopy();
				experiment.manual_primary_dist_hist.DrawCopy("same");
				ejech.DrawCopy("same");
			outfile->cd();
			can_view->SetName("CM_angles");can_view->Write("CM_angles",TObject::kOverwrite);
			gROOT->cd(); 			
			
			can_view->Clear();
			can_view->Divide(2);
			can_view->cd(1);
			experiment.mask_manual(t_d,t_u,true,false);//saves on wasted computation
				experiment.draw_hit_pattern_2D(1,500000,1,1,0);
				((TH2D*)gPad->GetPrimitive("Hit_Map_copy"))->SetTitle(r.c_str());
			can_view->cd(2);
				experiment.mask_manual(t_d,t_u,false,true);//saves on wasted computation
				experiment.draw_hit_pattern_2D(1,500000,1,0,1);
				((TH2D*)gPad->GetPrimitive("Hit_Map_copy"))->SetTitle(e.c_str());
			outfile->cd();
			can_view->SetName("Hit_map_ditribution");can_view->Write("Hit_map_ditribution",TObject::kOverwrite);
			gROOT->cd(); 		
			
			
			can_view->Clear();
			can_view->cd();
				experiment.mask_manual(t_d,t_u,true,true);//saves on wasted computation
				experiment.auto_E_theta(500000,1,1,1);
			outfile->cd();
			can_view->SetName("Theta_E");can_view->Write("Theta_E",TObject::kOverwrite);
			gROOT->cd(); 

		}
		
		if(control==5){spice_auto_setup(experiment,config,S3,overr);
			can_view->Clear();
			can_view->cd();
			experiment.draw_hits_3D(2);
		}
	}
	

	
	
// 
// 	// Next we spend 20 seconds drawing events so that if something is very wrong it can be corrected
// 	// We set multiplicity=1 so it will only show us S3 good events
// 	experiment.draw_hits_3D(2,1,true,0.5,1,true,20);
// 	//draw_hits_3D(projection,hit_multiplicity,obstructions,display_time,refresh_rate,draw_misses,run_time)
// 
// 	
// 	


	
	
	
	
if(save_text){
	std::cout.rdbuf(backup);        // restore cout's original streambuf
	filestr.close();
		
	string x;
	ifstream inFile(file_title+".txt");
	while (getline(inFile,x)){
	cout << x << endl ;
	}
	inFile.close();

	outfile->cd();
	TMacro mac((file_title+".txt").c_str()); 
	mac.Write("terminal.txt",TObject::kOverwrite);
}	


	outfile->Close();
	cout<<endl<<"/////////////////////////////////////////////////////////////////////////"<<endl;	
	app->Run();
	return 0;
}	
	
	
