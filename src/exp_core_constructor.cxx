
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


//////////////////////////////////////
/////////// Constructor   ///////////
/////////////////////////////////////

exp_core::exp_core(double worldsize_in) ://worldsize in mm
//////////////////////////////////////////////
/////////// Public Data Members    ///////////
//////////////////////////////////////////////
// C++11 CODE //part_names{"Recoil","Ejectile","Decay A","Decay B"},
// C++11 CODE //projection_names{"xyz","yz","yx","zx"},
dist_target_KE_beam(),dist_target_KE_0(),dist_target_beta(),dist_target_P_0(),dist_target_E_star(),dist_target_KE_1(),dist_target_P_1(),
density_targ(1.0),density_back(1.0),
implant_mode(1),
fusion_evaporation_simplify(0),
//////////////////////////////////////////////
/////////// Private Data Members    //////////
//////////////////////////////////////////////
current_target_fraction(0.5),decay_target_fraction(0.5),dist_target("dist_target","1"),dist_target_var(0),thick_target_calcs(0),stopper_seperation(0.0),
thick_target_fraction_0(0.5),thick_beta_CoM(0.0),thick_KE_0_tot_CoM(0.0),thick_P_0_CoM(0.0),thick_reco_E_star(0.0),thick_KE_1_tot_CoM(0.0),thick_P_1_CoM(0.0),
worldsize(worldsize_in),
targ_Z(6),targ_A(12),beam_Z(1),beam_A(1),targ_mass(12.0),beam_mass(1.0),targ(6,12,0),
ejec_Z(0),ejec_A(0),reco_Z(0),reco_A(0),ejec_mass(0.0),reco_mass(0),reco_E_star(0.0),
QVal(0.0),QMan(0.0),E_beam_min(0.0),basic_barrier_in(0.0),basic_barrier_out(0.0),E_beam_barrier(0.0),beta_CoM(0.0),
E_beam(0.0),P_beam(0.0),E_beam_targ(0.0),P_beam_targ(0.0),
KE_0_tot_CoM(0.0),P_0_CoM(0.0),KE_1_tot_CoM(0.0),P_1_CoM(0.0),
validQ(0),targ_fusion(0),
lorR(),lorE(),lorA(),lorB(),
offset_primary(0.0,0.0,0.0),offset_decay(0.0,0.0,0.0),offset_stopper(0.0,0.0,0.0),use_offsets(0),lifetime(0.0),post_beta(0.0),
simm_dist(0),simm_dist_parA(1),simm_dist_parB(-1),simm_dist_parC(0),reverse_primary(0),
decay_events(0),decay_type(0),
decay_Z_A(0),decay_Z_B(0),decay_A_A(0),decay_A_B(0),	
mass_decay1(0.0),mass_decay2(0.0),P_decay_recoil_CoM(0.0),KE_decay_tot_CoM(0.0),decay_recoil_E_star(0.0),
mass_ratio(0.0),delta_m(0.0),delta_TKE(0.0),ratioAZ(0.0),frag_bar_max(0.0),TKECentral(0),TKE_fiss_zero(0.0)
///////////////////////////////////
/////////// Function    ///////////
///////////////////////////////////
{
	rand.SetSeed();
	valid_doubles.resize(0);
	detectors.resize(0);
	
	valid_dets.assign(4,vector< bool >());
	det_hits.assign(4,vector< int >());
	stop_distance.assign(4,-1);
	
	offsets.resize(0);
	offsets.push_back(&offset_primary);
	offsets.push_back(&offset_primary);
	offsets.push_back(&offset_decay);
	offsets.push_back(&offset_decay);
	
	lor_point.resize(0);
	lor_point.push_back(&lorR);
	lor_point.push_back(&lorE);
	lor_point.push_back(&lorA);
	lor_point.push_back(&lorB);
	
	part_Z.resize(0);
	part_Z.push_back(&reco_Z);
	part_Z.push_back(&ejec_Z);
	part_Z.push_back(&decay_Z_A);
	part_Z.push_back(&decay_Z_B);
	
	part_A.resize(0);
	part_A.push_back(&reco_A);
	part_A.push_back(&ejec_A);
	part_A.push_back(&decay_A_A);
	part_A.push_back(&decay_A_B);
	
	part_M.resize(0);
	part_M.push_back(&reco_mass);
	part_M.push_back(&ejec_mass);
	part_M.push_back(&mass_decay1);
	part_M.push_back(&mass_decay2);

	//initilise some graphics
	threedee= new TH3D("threedee","threedee",30,-worldsize,worldsize,30,-worldsize,worldsize,30,-worldsize,worldsize);
	threedee->GetXaxis()->SetTitle("X-axis [mm]");
	threedee->GetYaxis()->SetTitle("Z-axis [mm]");
	threedee->GetZaxis()->SetTitle("Y-axis [mm]");
	
	twodee= new TH2D("twodee","twodee",100,-worldsize,worldsize,100,-worldsize,worldsize);
	
	tp_proj= new TH2D("tp_proj","tp_proj",360,0,pi,720,0,2*pi);
	tp_proj->GetXaxis()->SetTitle("#Theta_{LAB}");
	tp_proj->GetYaxis()->SetTitle("#Phi_{LAB}");	
	
	for(int i=0;i<4;i++){
		part_pl3.push_back(new TPolyLine3D(2));
		part_pl.push_back(new TPolyLine(2));
		part_pl3[i]->SetLineColor(i+2);
		part_pl3[i]->SetLineWidth(2);
		part_pl[i]->SetLineColor(i+2);
		part_pl[i]->SetLineWidth(2);
	}

	//initilise thetaM,phiM;
	this->gen_uniform(thetaM,phiM);
	this->gen_uniform(thetaS,phiS);
	
	dist_target_hist = TH1D("targ_interaction","targ_interaction",1000,0,1);
	beta_dist_hist = TH1D("beta_dist_hist","beta_dist_hist",1000,0,1);
	manual_primary_dist_hist= TH1D("manual_primary_dist_hist","manual_primary_dist_hist",5000,0,pi);
	manual_primary_dist_hist.Copy(manual_primary_dist_store);
	manual_primary_dist_hist.Copy(manual_primary_dist_mask);
	this->reset_mask();

	part_names.push_back("Recoil");
	part_names.push_back("Ejectile");
	part_names.push_back("Decay A");
	part_names.push_back("Decay B");

	projection_names.push_back("xyz");
	projection_names.push_back("yz");
	projection_names.push_back("yx");
	projection_names.push_back("zx B");

// C++11 CODE //{"","","",""},

}


//////////////////////////////////////
///////////  Destructor   ///////////
/////////////////////////////////////

exp_core::~exp_core()
{
	delete threedee;
	delete tp_proj;
	delete twodee;
	
	for(int i=0;i<4;i++){	
		delete part_pl[i];
		delete part_pl3[i];
	}
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
