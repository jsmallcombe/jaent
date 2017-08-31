
//
//
// James's not Geant4 : jaent V1.3.1
// james.smallcombe@outlook.com 01/4/2015
//
//


#include "exp_core.h"

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// Non-class free functions  ///////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// Public Members Functions  ///////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

//set basic particles

void exp_core::set_beam(int Z,int A,double MeV,double mass)
{
	beam_Z=Z;
	beam_A=A;
	if(mass>0) beam_mass=mass;
	else beam_mass=nuclear_data_ob::get_mass(Z,A);
	this->set_ejectile();	
	this->check_fusion();
	this->set_all_base_E_calcs();
	if(MeV>0) this->set_E_beam(MeV);
}
void exp_core::set_beam(string S,int A,double MeV,double mass)
{this->set_beam(nuclear_data_ob::get_Z(S),A,MeV,mass);}



void exp_core::set_targ(target T,double mass)
{	thick_target_calcs=false;
	stopper_seperation=0.0;
	targ=T;
	targ_Z=T.Z();
	targ_A=T.A();
	density_targ=nuclear_data_ob::get_nom_density(T.Z());
	density_back=nuclear_data_ob::get_nom_density(T.bZ());
	if(mass>0) targ_mass=mass;
	else targ_mass=nuclear_data_ob::get_mass(targ_Z,targ_A);	
	this->set_ejectile();	
	this->check_fusion();	
	this->set_all_base_E_calcs();
}
void exp_core::set_targ(int Z,int A,double mass){this->set_targ(target(Z,A,0),mass);}
void exp_core::set_targ(string S,int A,double mass){this->set_targ(nuclear_data_ob::get_Z(S),A,mass);}


void exp_core::set_reco(int Z,int A,double mass)
{
	reco_Z=Z;
	reco_A=A;
	if(mass>0) reco_mass=mass;
	else reco_mass=nuclear_data_ob::get_mass(Z,A);	
	this->set_ejectile();	
	this->check_fusion();
	this->set_all_base_E_calcs();
}
void exp_core::set_reco(string S,int A,double mass)
{this->set_reco(nuclear_data_ob::get_Z(S),A,mass);}


void exp_core::set_ejec(int Z,int A,double mass)
{
	ejec_Z=Z;
	ejec_A=A;
	if(mass>0) ejec_mass=mass;
	else ejec_mass=nuclear_data_ob::get_mass(Z,A);	
	
	//set_ejectile
	reco_Z=targ_Z+beam_Z-ejec_Z;
	reco_A=targ_A+beam_A-ejec_A;
	reco_mass=nuclear_data_ob::get_mass(reco_Z,reco_A);
		
	this->check_fusion();
	this->set_all_base_E_calcs();
}
void exp_core::set_ejec(string S,int A,double mass)
{this->set_ejec(nuclear_data_ob::get_Z(S),A,mass);}


void exp_core::set_reaction(int Z1,int A1,int Z2,int A2,int Z3,int A3,double beE)
{
	this->set_targ(Z1,A1);
	this->set_beam(Z2,A2,beE);
	this->set_reco(Z3,A3);
	if(beE>0)this->set_E_beam(beE);
}
void exp_core::set_reaction(string S1,int A1,string S2,int A2,string S3,int A3,double beE)
{this->set_reaction(nuclear_data_ob::get_Z(S1),A1,nuclear_data_ob::get_Z(S2),A2,nuclear_data_ob::get_Z(S3),A3,beE);}

/// set energy modification

void exp_core::set_E_beam(double E)
{
	E_beam=E;
	this->set_beam_calcs();
}

void exp_core::set_E_cm_beam(double E){
	double Eb= reverseal_mom_calc_KE(E,beam_mass,targ_mass);
	Eb= get_KE(Eb,beam_mass);
	if(targ.GetThickness()>0){
		if(!thick_target_calcs) Eb= targ.beam_e_centre_reverse(beam_Z,beam_A,Eb);
		Eb= targ.beam_e_centre_reverse(beam_Z,beam_A,Eb,TVector3(0.0,0.0,1.0),thick_target_fraction_0);
	}
	this->set_E_beam(Eb);
}


void exp_core::set_fusion(double mass){
	targ_fusion=true;	
	reco_Z=targ_Z+beam_Z;
	reco_A=targ_A+beam_A;
	if(mass>0) reco_mass=mass;
	else reco_mass=nuclear_data_ob::get_mass(reco_Z,reco_A);	
	this->set_ejectile();	
	this->set_QV();	
	this->set_beam_calcs();
	this->set_calc_threshold();	
}

void exp_core::set_elastic(){
	this->set_Q_manual(0.0);
	reco_Z=targ_Z;
	reco_A=targ_A;
	reco_mass=targ_mass;
	ejec_Z=beam_Z;
	ejec_A=beam_A;
	ejec_mass=beam_mass;
	reco_E_star=0;
	this->set_all_base_E_calcs();
	targ_fusion=false;
}	

void exp_core::set_Q_manual(double Qma)
{
	QMan=Qma;
	QVal=Qma;
	validQ=true;
	this->set_beam_calcs();
	this->set_calc_threshold();	
}

void exp_core::set_E_star(double E)
{
	reco_E_star=abs(E);
	this->set_calc_threshold();
	if(targ_fusion)E_beam=E_beam_min;	
	this->set_beam_calcs();
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void exp_core::set_all_base_E_calcs(){
	this->set_QV();
	this->set_beam_calcs();
	this->set_barrier();
	this->set_calc_threshold();
	fusion_evaporation_simplify=false;
}


void exp_core::set_ejectile(){
	ejec_Z=targ_Z+beam_Z-reco_Z;
	ejec_A=targ_A+beam_A-reco_A;
	if(ejec_A<1){ejec_A=0;ejec_Z=0;}
	ejec_mass=nuclear_data_ob::get_mass(ejec_Z,ejec_A);
	
}

void exp_core::set_QV(){
	if(QMan!=0.0){
		QVal=QMan;
	}else{
		QVal=(beam_mass+targ_mass-ejec_mass-reco_mass)*jam_phys_amu;	
		QVal=round(QVal*1000)/1000.0;//strange rounding error before, so now to nearest KeV

		if(nuclear_data_ob::is_mass(beam_Z,beam_A)&& nuclear_data_ob::is_mass(targ_Z,targ_A)
		&& nuclear_data_ob::is_mass(ejec_Z,ejec_A)&& nuclear_data_ob::is_mass(targ_Z+beam_Z-ejec_Z,beam_A+targ_A-ejec_A)){
			validQ=true;	
		}else{
			validQ=false;
		}
	}

}

void exp_core::set_beam_calcs(){
	current_target_fraction=0.5;decay_target_fraction=0.5;
	
	if(targ.GetThickness()>0)E_beam_targ=targ.beam_e_centre(beam_Z,beam_A,E_beam);
	else E_beam_targ=E_beam;
	P_beam_targ=get_rel_mom(E_beam_targ,beam_mass);

	beta_CoM=get_com_beta(P_beam_targ,beam_mass,targ_mass);
	KE_0_tot_CoM=get_com_KE(P_beam_targ,beam_mass,targ_mass);	
	P_0_CoM=momentum_beam_P_CoM(P_beam_targ,beam_mass,targ_mass);
	
	if(targ_fusion){
		reco_E_star=KE_0_tot_CoM+QVal;
		if(reco_E_star<0)reco_E_star=0;		
	}	
	KE_1_tot_CoM=KE_0_tot_CoM+QVal-reco_E_star;
	if(KE_1_tot_CoM<0)KE_1_tot_CoM=0;
	
	P_1_CoM=momentum_energysplit_CoM(KE_1_tot_CoM,(reco_mass+(reco_E_star/jam_phys_amu)),ejec_mass);
	
	if(targ.GetThickness()>0&&thick_target_calcs)this->calc_thick_target_interaction();
}

void exp_core::set_barrier(){
	basic_barrier_in=classical_barrier(targ_A,targ_Z,beam_A,beam_Z);
	basic_barrier_out=classical_barrier(ejec_A,ejec_Z,reco_A,reco_Z);
	
	double P_min=reverseal_mom_calc_KE(basic_barrier_in,beam_mass,targ_mass);
	E_beam_barrier=get_KE(P_min,beam_mass);
// 	if(targ.GetThickness()>0)E_beam_barrier=targ.beam_e_centre_reverse(beam_Z,beam_A,E_beam_barrier,TVector3(0.0,0.0,1.0),current_target_fraction);

	if(targ.GetThickness()>0){
		if(!thick_target_calcs) E_beam_barrier=targ.beam_e_centre_reverse(beam_Z,beam_A,E_beam_barrier);
		else E_beam_barrier=targ.beam_e_centre_reverse(beam_Z,beam_A,E_beam_barrier,TVector3(0.0,0.0,1.0),thick_target_fraction_0);
	}
}

void exp_core::set_calc_threshold(){
	double Qtot=QVal-reco_E_star;
	if(Qtot>0)Qtot=0;
	double P_min=reverseal_mom_calc_KE(abs(Qtot),beam_mass,targ_mass);
	E_beam_min=get_KE(P_min,beam_mass);
	if(targ.GetThickness()>0){
		if(!thick_target_calcs) E_beam_min=targ.beam_e_centre_reverse(beam_Z,beam_A,E_beam_min);
		else E_beam_min=targ.beam_e_centre_reverse(beam_Z,beam_A,E_beam_min,TVector3(0.0,0.0,1.0),thick_target_fraction_0);
	}
}

void exp_core::check_fusion(){
	if(((targ_Z+beam_Z==reco_Z)&&(targ_A+beam_A==reco_A))||
		((targ_Z+beam_Z==ejec_Z)&&(targ_A+beam_A==ejec_A)))
		this->set_fusion();
	else	{targ_fusion=false;reco_E_star=0;}
}
