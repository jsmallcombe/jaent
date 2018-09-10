
//
//
// James's not Geant4 : jaent V1.3.1
// james.smallcombe@outlook.com 01/4/2015
//
//




#ifndef __EXPCORE_H_INCLUDED__   // if x.h hasn't been included yet...
#define __EXPCORE_H_INCLUDED__   //   #define this so the compiler knows it has been included

#include <TSystem.h>


#include <iostream>
#include <iomanip> 
#include <fstream>
#include <vector>
#include <stdlib.h>     //for using the function "sleep"
#include <sstream>

// #include <algorithm>
// #include <string>
// #include <sstream>
// #include <queue>
// #include <list>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTree.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TTreePlayer.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRotation.h>
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>

#include <TCutG.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TRandom2.h>
#include <TText.h>
#include <TStopwatch.h>
#include <TFormula.h>
#include <TVirtualPad.h>
#include <TPad.h>
#include <TStyle.h>

#include <jlibmaster.h>
#include <jphysmaster.h>
// #include "james_nuclear_data_ob.h"
// #include "james_detector_func.h"
// #include "james_physics.h"//includes root_maths
// #include "james_fission.h"
// #include "james_target.h"


#include "detector_class.h"


using namespace std;


// Should always call set particles BEFORE manually setting Q value

/////////////////////////////////////////////////////////////////////////////////////
/////////// "exp_core" The core class for holding experimental setups     ///////////
///////////           and performinge experiment planning simms           ///////////  
/////////////////////////////////////////////////////////////////////////////////////
class exp_core
{

//////////////////////////////////////////
//////////////////////////////////////////
///////////     Data Members   ///////////
//////////////////////////////////////////
//////////////////////////////////////////
public:
	vector< string > part_names;
	vector< string > projection_names;
	//thick_target stuff
	TH1D dist_target_hist,beta_dist_hist,manual_primary_dist_hist,manual_primary_dist_store,manual_primary_dist_mask;
	TGraph dist_target_KE_beam,dist_target_KE_0,dist_target_beta,dist_target_P_0,dist_target_E_star,dist_target_KE_1,dist_target_P_1;
	double density_targ,density_back;
	int implant_mode;

	bool fusion_evaporation_simplify;
	
	int get_TZ(){return targ_Z;};
	int get_TA(){return targ_A;};
	int get_BZ(){return beam_Z;};
	int get_BA(){return beam_A;};
	int get_RZ(){return reco_Z;};
	int get_RA(){return reco_A;};
	int get_EZ(){return ejec_Z;};
	int get_EA(){return ejec_A;};
	int get_D1Z(){return decay_Z_A;};
	int get_D1A(){return decay_A_A;};
	int get_D2Z(){return decay_Z_B;};
	int get_D2A(){return decay_A_B;};
	
	//Useful little access function if doing a manual loop externally
	TLorentzVector GetLor(int i){if(abs(i)<4)return *lor_point[i];else return TLorentzVector();}
	vector< vector< int > > get_det_hits(){return det_hits;}
	
private:	
	//thick_target stuff
	double current_target_fraction,decay_target_fraction;
	TFormula dist_target;
	int dist_target_var;
	bool thick_target_calcs; 
	double stopper_seperation;
	double thick_target_fraction_0,thick_beta_CoM,thick_KE_0_tot_CoM,thick_P_0_CoM,thick_reco_E_star,thick_KE_1_tot_CoM,thick_P_1_CoM;
	
	//graphics
	TH3D	*threedee;
	TH2D	*tp_proj,*twodee; 
	vector<  TPolyLine3D* > part_pl3;
	vector<  TPolyLine* > part_pl;
	double worldsize;
	    
	//beam target info
	int targ_Z,targ_A,beam_Z,beam_A;
	double targ_mass,beam_mass;
	target targ;

	//product info
	int ejec_Z,ejec_A,reco_Z,reco_A;
	double ejec_mass;
	double reco_mass; //does not include excitation mass (so add where needed)
	double reco_E_star;

	//Energy data
	double QVal,QMan;
	double E_beam_min;
	double basic_barrier_in,basic_barrier_out;
	double E_beam_barrier;
	double beta_CoM;
	double E_beam, P_beam;
	double E_beam_targ,P_beam_targ;
	double KE_0_tot_CoM,P_0_CoM,KE_1_tot_CoM,P_1_CoM;
	    
	//detectors
	vector< vector<bool> > valid_doubles; //if one hit hits two detectors, allowed?
	vector< detector > detectors;
	//detectors_hit_data
	//0 = Recoil, 1 = Ejectile, 2 = decay_A (recoil), 3 = decay_B (eject)
	vector< vector< bool > > valid_dets;
	vector< vector< int > > det_hits;
	vector< double > stop_distance;//doesnt work for gammas

	//other
	TRandom2 rand;
	bool validQ;
	bool targ_fusion;
	
	//random trajectory generation master variabls
	double thetaM,phiM,thetaS,phiS;
	
	//FINAL DATA LORENTZ VECTORS
	TLorentzVector lorR,lorE,lorA,lorB;
	vector< TLorentzVector* > lor_point;
	vector< int* > part_Z;
	vector< int* > part_A;
	vector< double* > part_M;
	vector< TVector3* > offsets;
	TVector3 offset_primary,offset_decay,offset_stopper;
	bool use_offsets;
	double lifetime;//g/cm^3
	double post_beta;

	
	//Control the CM distribution of the primary reaction
	int simm_dist;
	double simm_dist_parA,simm_dist_parB,simm_dist_parC;
	bool reverse_primary;
	
	//DECAY VARIALBES
	bool decay_events;
	int decay_type;
	int decay_Z_A,decay_Z_B,decay_A_A,decay_A_B;	
	double mass_decay1,mass_decay2;	
	double P_decay_recoil_CoM,KE_decay_tot_CoM;
	double decay_recoil_E_star;
	
		//Fission 
		double mass_ratio;
		double delta_m;
		double delta_TKE;
		double ratioAZ;
		double frag_bar_max;
		bool TKECentral;
		double TKE_fiss_zero;
		
		//Gamma
		
		//b_ray
		
		//ICE
		
		//alpha
	
	
////////////////////////////////////
////////////////////////////////////
//////////     Methods   ///////////
////////////////////////////////////
////////////////////////////////////	
	
/////////////////////////////////////////////////
/////////// Constructor & Destructor  ///////////
/////////// exp_core_constructor.cxx  ///////////
/////////////////////////////////////////////////	
public:
	// constructor
	exp_core(double=400);
	// destructor
	~exp_core();
	
///////////////////////////////////////////////
///////////   Detector management   ///////////
/////////// exp_core_detectors.cxx  ///////////
///////////////////////////////////////////////
public:
	void add_detector(vector< vector<double> >,TVector3=TVector3(0.0,0.0,1.0),TRotation=TRotation(),bool=false,bool=false,bool=false,bool=false);
	void add_detector(vector<double>,vector<double>,TVector3=TVector3(0.0,0.0,1.0),TRotation=TRotation(),bool=false,bool=false,bool=false,bool=false);
	void add_detector(detector,bool=false,bool=false,bool=false,bool=false);	
	
	void set_pair(int,int,bool=true);
	void reset_detectors();
	void reset_valid_dets();
	void set_implant_escape(int in=-1){implant_mode=in;}
	//0=no implants 1(initial)=not recorded in implant -1=recorded but needs detector pairing
	
	void set_valid_dets(int,int);
	void set_valid_part(int in,bool=true,bool=true,bool=true,bool=true);

	
	void reset_doubles();	
	int detN(){return detectors.size();}
	void SetDetName(string name){SetDetName(name,detN()-1);}
	void SetDetName(string name,unsigned int d){if(d<(unsigned)detN())detectors[d].name=name;}
	string GetDetName(unsigned int d){
		if(d<(unsigned)detN()){
			stringstream ss;
			ss<<detectors[d].name;
			if(ss.str().size()<1)ss<<"Detector "<<d;
			return ss.str();
		}
		return "";
	}
	
	//intentionally a copy rather than a pointer because detectors are sacred
	detector get_det(unsigned int in=0){if(in<(unsigned)this->detN())return detectors[in]; else return detector();}
private:
	bool det_hit_quick(int,int);
	void det_check_hits(int,bool);
	void det_check_hits_all(bool);
	
	bool all_hit_quick(double theta ,double phi);
	
	int current_multiplicity();
	
	double max_distance_hit(int);
	
	void addhitsdet();	//add all hits to their detectors run after check hits
	void extend_doubles();
//////////////////////////////////////////////
///////////  Miscellaneous Fn    ///////////
///////////   exp_core_misc.cxx    ///////////
//////////////////////////////////////////////
public:
	void set_gun(int=6,int=12,double=20.0);
	void set_gun(int,int,double,target);
	void spontaniously_decay();
	
	//get current physical parameters

	double get_reco_E_star(){this->null_target_pos(); return reco_E_star;}
	double get_E_beam_min(){this->null_target_pos(); return E_beam_min;}
	double get_basic_barrier_in(){this->null_target_pos(); return basic_barrier_in;}
	double get_E_beam_barrier(){this->null_target_pos(); return E_beam_barrier;}
	double get_beta_CoM(){this->null_target_pos(); return beta_CoM;}
	double get_E_beam(){this->null_target_pos(); return E_beam;}
	double get_E_beam_targ(){this->null_target_pos(); return E_beam_targ;}
	double get_KE_0_tot_CoM(){this->null_target_pos(); return KE_0_tot_CoM;}
	double get_KE_1_tot_CoM(){this->null_target_pos(); return KE_1_tot_CoM;}
	double get_P_1_CoM(){this->null_target_pos(); return P_1_CoM;}

private:		
/////////////////////////////////////////////////
/////////// Screen Printing check setup /////////
///////////    exp_core_printing.cxx    /////////
/////////////////////////////////////////////////
public:
	void print_reaction();
	void print_decay();
	void print_detectables();
	void print_doubles();
	void print_target();
	string channelnames(int,int);
private:
////////////////////////////////////////////////////////
///////////   Set core reaction parameters   ///////////
/////////// 	  exp_core_par_set.cxx       ///////////
////////////////////////////////////////////////////////
public:
	void set_beam(int,int,double=0.0,double=0.0);
	void set_targ(int,int,double=0.0);
	void set_targ(target,double=0.0);
	void set_reco(int,int,double=0.0);
	void set_ejec(int,int,double=0.0);
	
	//stringinput
	void set_beam(string,int,double=0.0,double=0.0);
	void set_targ(string,int,double=0.0);
	void set_reco(string,int,double=0.0);
	void set_ejec(string,int,double=0.0);	
	
	void set_reaction(int,int,int,int,int,int,double=0.0);	
	void set_reaction(string,int,string,int,string,int,double=0.0);
	
	void set_E_beam(double);
	void set_E_cm_beam(double);
	void set_fusion(double=0);	
	void set_elastic();	
	void set_Q_manual(double);	
	void set_E_star(double);
private:		
	void set_all_base_E_calcs();
	void set_ejectile();
	void set_QV();
	void set_beam_calcs();
	void set_barrier();
	void set_calc_threshold();		
	void check_fusion();	
	
//////////////////////////////////////////////////////
///////////     Main event generation      ///////////
///////////     exp_core_event_simm.cxx    ///////////
//////////////////////////////////////////////////////
public:
	void basic_make_single(bool =false,int=0,int=0,int=0,int=0);
	void basic_hit_count(int,bool =false,int=0,int=0,int=0,int=0);
	double basic_det_hit_frac(int,int,bool=false);
	double basic_det_hit_multi(int,int,bool=false);
	void basic_decay_check(int);
private:		

	void gen_event();
	void gen_primary();
	void gen_decay();
	void gen_Eloss();

//////////////////////////////////////////////////////
///////////     Primary Distributions      ///////////
///////////     Rutherford Scattering      ///////////
///////////   exp_core_distributions.cxx   ///////////
//////////////////////////////////////////////////////
public:
	
	
	double set_uniform(double=0,double=pi);//return fraction 
	void set_zeroid(double=0.1);
	void set_primary_dist(string);
	void set_primary_dist(TGraph*);
	void set_primary_dist(TFormula*);
	void set_primary_dist(TF1*);
	void set_ruthmask(double,double,bool=true,bool=true);
	double* mask_manual(double,double,bool=true,bool=true);
	
	double frac_masked(){return manual_primary_dist_hist.Integral();}//Ask what fraction we have masked too, as unmask integral = 1	
	
	double set_rutherford(double=0.1,double=pi);
	double set_rutherford_lab(double=0.1,double=pi);
	// Set simm to rutherford scattering over given range
	// Returns integrated phi2pi rutherford_crosssection in mb [Inputs(thetamin_cm,thetamax_cm)]
	double get_rutherford_lab(double=0.1,double=pi,bool=true,bool=true);
	
	void reverse_primary_dist(){reverse_primary=!reverse_primary;this->apply_mask();}
	
	void auto_rutherford(int=10000,bool=false);
	void auto_rutherford_OTT(int=10000,bool=false);//do backing aswell
private:	
	void gen_uniform(double& ,double&);
	void gen_ruth(double& ,double&);
	void gen_manual_primary(double& ,double&);
	void gen_zeroid(double& ,double&);
	void reset_mask();
	void apply_mask();
	double* generate_mask(double,double,bool,bool);
	
	// 0 uniform: min angle, max angles, returns 2pi fraction
	// 1 rutherford: min angle, max angles, returns mb crosssection
	// 2 zeroid, gausian of theta_cm centred at 0, sigma in radians input
	
////////////////////////////////////////////////////////////
//////////       Decay reaction functions        ///////////
//////////   Parameters, event code and outputs  ///////////
//////////          exp_core_decay.cxx           ///////////
////////////////////////////////////////////////////////////
public:
	void set_gamma(double =0.0);
	void set_b_ray(bool=false,double=0.0,double =0.0);	
	void set_ICE(double=0.0);
	void set_alpha(double=0.0,double=0.0);
	void decay_off(){decay_events=false;}
private:
////////////////////////////////////////////////////////////
//////////       Decay fission functions         ///////////
//////////   Parameters, event code and outputs  ///////////
//////////          exp_core_fission.cxx         ///////////
////////////////////////////////////////////////////////////
public:
	void set_fission(double=1.0,double=0.0,double=0.0,int=0);
	TH2D* fission_fragment_dist();
	TH2D* fission_fragment_dist_detect(bool=false);
	TH2D* fission_fragment_angles(int=500000);
	TH2D* fission_fragment_angles_detect(bool=false,int=500000);
	TH3D* fission_fragment_transfer_angles(bool=false);
	TH1D* fission_fragment_dT(int detA,int detB,bool=false);
private:
	void gen_fission();

///////////////////////////////////////////////////////
///////////     Graphical output functions  ///////////
///////////         exp_coredrawing.cxx         ///////////
///////////////////////////////////////////////////////
public:
	//drawing
	void draw_exp(int=0,int=1);//best already be pointing to a pad
	void draw_phi(bool=false);//best already be pointing to a pad
	void draw_xz_labels();//best already be pointing to a pad
	void draw_boost_detectors(bool=true);//best already be pointing to a pad
	void draw_target_interaction(int mult=0,bool obstruc=false,int reps=10000000);//best already be pointing to a pad
	void draw_decay_Z(int=0,bool=false,double=0);//best already be pointing to a pad
	
	
	TH1D theta_cover_auto_norm(bool=false,int=1000000);
	TH1D theta_cover_auto(bool=false,int=1000000);
	TH1D theta_hist(int=1000000,bool=false,bool=false,bool=true,bool=false,bool=false);
	
	void draw_primary_kinematics();//best already be pointing to a pad
	
	void draw_hit_pattern_2D(int=1,int=10000,bool=1,bool=1,bool=0,bool=0,bool=0);//best already be pointing to a pad
	void auto_E_theta(int=10000,bool=1,bool=1,bool=0,bool=0,bool=0);//best already be pointing to a pad

	void draw_hits_3D(int=0,int=1,bool=true,double=0.75,int=1,bool=true,double=120);
//void draw_hits_3D(projection,particle_hit_multiplicity,obstructions,display_time,refresh_rate,counterpart_misses,run_time)

	
private:	
	void draw_arrow();
	TVector3 path_vec_calc(int,int det=-1);	
	void draw_path_3D(int,int=0,int=-1);
	void draw_fill(int=0);
	void kill_lines();
	
/////////////////////////////////////////////////////////////////////////
///////////     Calculations for thick target interactions    ///////////
///////////            exp_core_thick_target.cxx              ///////////
/////////////////////////////////////////////////////////////////////////
public:
	void set_target_interaction_off();
	void set_target_interaction(int=0);
	void set_target_interaction(TFormula,int=0);
private:
	void calc_thick_target_interaction();
	void gen_random_target_pos();
	void null_target_pos();
	void lab_loretz_targ_exit(int,double);
	
/////////////////////////////////////////////////////////////////////////
///////////     Calculations for offsets and lifetimes        ///////////
///////////           exp_core_offset_lifetime.cxx            ///////////
/////////////////////////////////////////////////////////////////////////
public:
	void set_target_primary_offset(TVector3=TVector3(0.0,0.0,0.0));
	void set_target_primary_offset(double,double,double);
	void set_lifetime_ns(double=0.0,double=0.0,double=0.0);
	void set_halflife_ns(double=0.0,double=0.0,double=0.0);
	void set_stopper_seperation(double sep){stopper_seperation=abs(sep);}
private:
	
	double get_range(int,bool,double);
	double passage(int,bool,double,double);
	double lifetime_track(int ,double& ,double&);

	
	
	friend void add_target_ladder(exp_core&,double,double,double);//needs to be friend to fetch target_norm
	
	
////////////////////////////////////////////////////////
///////////  gosia calculations and loops    ///////////
////////////////////////////////////////////////////////
public:
	void BuildOPINTI(int Nmesh=20,int NIntegrate=60,double thetamax=0);	
	void ReadDrawMesh(string filename="OPINTI.txt",double Z=32.);
private:
	void FindThetaMinMax(double& tMin,double& tMax);
	void ShapeToPhiList(vector<double>& ThetaList,vector<vector<double>>& PhiList,int NPoints,double TMin,double TMax);
	
};










#endif 
