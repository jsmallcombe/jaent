
//
//
// James's not Geant4 : jaent V1.3.1
// james.smallcombe@outlook.com 01/4/2015
//
//


#include "exp_core.h"

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// Public Members Functions  ///////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////


void exp_core::set_gamma(double keV){
	double MeV=keV/1000.0;
	if(MeV<=0)MeV=reco_E_star;
	
	decay_events=true;
	decay_type=1;
	
	mass_decay1=reco_mass;
	mass_decay2=0.0;	
	decay_A_A=reco_A;
	decay_A_B=0;
	decay_Z_A=reco_Z;
	decay_Z_B=0;
	
	if(reco_E_star<MeV)reco_E_star=MeV;
	
	P_decay_recoil_CoM=momentum_energysplit_CoM(MeV,mass_decay1,mass_decay2);
	decay_recoil_E_star=reco_E_star-MeV;
}



void exp_core::set_b_ray(bool betaplus,double T_MeV,double final_MeV){	
	
	decay_events=true;
	decay_type=2;
	
	int plus=-1;
	if(betaplus)plus=1;

	decay_A_A=reco_A;
	decay_Z_A=reco_Z-plus;	
	decay_A_B=0;
	decay_Z_B=plus;	
	mass_decay2=0.511/jam_phys_amu;	
	
	if(final_MeV>0){
		mass_decay1=reco_mass-(0.511+T_MeV+final_MeV-reco_E_star)/jam_phys_amu;
		decay_recoil_E_star=final_MeV;
	}else{
		mass_decay1=nuclear_data_ob::get_mass(decay_Z_A,decay_A_A);
		double Q=(reco_mass-mass_decay1)*jam_phys_amu-0.511+reco_E_star;//total available
		
		if(Q>T_MeV&&T_MeV>0){
			decay_recoil_E_star=Q-T_MeV;//going to excited state
		}else{
			if(T_MeV<=0)T_MeV=0.001;
			//overrides the mass
			decay_recoil_E_star=0;
			mass_decay1=reco_mass-(0.511+T_MeV-reco_E_star)/jam_phys_amu;
		}
	}
	
	beta_dist_hist=TH1D("beta_dist_hist","beta_dist_hist",1000,0,T_MeV);
	for(int i=1;i<=beta_dist_hist.GetNbinsX();i++){
		double x=beta_dist_hist.GetXaxis()->GetBinCenter(i);
		beta_dist_hist.SetBinContent(i,unnorm_beta_dist_TKE(x,T_MeV,decay_Z_A*(-plus)));
	}
	beta_dist_hist.Scale(1/beta_dist_hist.Integral());
	beta_dist_hist.ComputeIntegral();
	
}

void exp_core::set_ICE(double keV){	//basically gamma
	decay_events=true;
	decay_type=3;

	double MeV=keV/1000;
	double bind=K_bind_aprox_keV(reco_Z)/1000.0;
	if(MeV<=0)MeV=reco_E_star;

	if(MeV<=bind){bind=0;cout<<endl<<"ICE transition energy less than binding."<<endl;}
	
	mass_decay1=reco_mass;
	mass_decay2=0.511/jam_phys_amu;	
	decay_A_A=reco_A;
	decay_A_B=0;
	decay_Z_A=reco_Z;
	decay_Z_B=-1;
	
	if(reco_E_star<MeV)reco_E_star=MeV;
	
	P_decay_recoil_CoM=momentum_energysplit_CoM(MeV-bind,mass_decay1,mass_decay2);
	decay_recoil_E_star=reco_E_star-MeV;
}



void exp_core::set_alpha(double T_MeV,double final_MeV){ // basically beta
	decay_events=true;
	decay_type=4;

	decay_A_A=reco_A-4;
	decay_Z_A=reco_Z-2;	
	decay_A_B=4;
	decay_Z_B=2;	
	mass_decay2=nuclear_data_ob::get_mass(2,4);	
	
	if(final_MeV>0){
		mass_decay1=reco_mass-mass_decay2-(T_MeV+final_MeV-reco_E_star)/jam_phys_amu;
		decay_recoil_E_star=final_MeV;
	}else{
		mass_decay1=nuclear_data_ob::get_mass(decay_Z_A,decay_A_A);
		double Q=(reco_mass-mass_decay1-mass_decay2)*jam_phys_amu+reco_E_star;//total available
		
		if(Q>T_MeV&&T_MeV>0){
			decay_recoil_E_star=Q-T_MeV;//going to excited state
		}else{
			if(T_MeV<=0)T_MeV=0.001;
			decay_recoil_E_star=0;
			mass_decay1=reco_mass-mass_decay2-(T_MeV-reco_E_star)/jam_phys_amu;
		}
	}
	
	P_decay_recoil_CoM=momentum_energysplit_CoM(T_MeV,mass_decay1,mass_decay2);
}



/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

// void exp_core::(){
// }
// void exp_core::gen_gamma(){
// 		P_decay_recoil_CoM=momentum_energysplit_CoM(KE_decay_tot_CoM,mass_decay1,mass_decay2);
// }
