
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

void exp_core::set_target_interaction_off()
{thick_target_calcs=false;this->set_beam_calcs();}

//some default cases
void exp_core::set_target_interaction(int select){
	TFormula distribution;
	int var=0;
	switch ( select ){
		case 1:	distribution=TFormula("dist_target","x");var=0;
			break;
		case 2:	distribution=TFormula("dist_target","x*x");var=0;
			break;
		case 3:	distribution=TFormula("dist_target","1/(x*x)");var=0;
			break;
		case 4:	distribution=TFormula("dist_target","x*(x>[0])+exp(x-[0])*x*(x<=[0])");distribution.SetParameter(0,basic_barrier_in);var=0;
			break;
		case 5:	distribution=TFormula("dist_target","(x>0)*[1]*[1]/(((x-[0])*(x-[0]))+[1]*[1])");distribution.SetParameters(4.5,0.1);var=3;
			break;
		default:distribution=TFormula("dist_target","1");
	}
	this->set_target_interaction(distribution,var);
}

// Interaction only happens in target material not backing
// Formula must be 1D and depend on (x) where var sets
// 0 -> x=KE_0_tot_CoM
// 1 -> x=P_0_CoM
// 2 -> x=beta_CoM
// 3 -> x=E_star
void exp_core::set_target_interaction(TFormula distribution,int var){
	dist_target_var=var;
	dist_target=distribution;
	thick_target_calcs=true;
	this->set_beam_calcs();
	this->set_barrier();
	this->set_calc_threshold();
}
		
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
	
	
//basically does "beam calcs" (from set_par) for points through the target
void exp_core::calc_thick_target_interaction(){

// // 	dist_target.Compile(); //apparently I didnt need this line, also it forgets parameters after so have to evaluate with TFormula::EvalPar 
	dist_target_hist.Reset();
	dist_target_KE_beam = TGraph();
	dist_target_KE_0 = TGraph();
	dist_target_beta = TGraph();		
	dist_target_P_0 = TGraph();
	dist_target_E_star = TGraph();
	dist_target_KE_1 = TGraph();
	dist_target_P_1 = TGraph();
	
	for(int i=1;i<=dist_target_hist.GetNbinsX();i++){
		double fraction=dist_target_hist.GetXaxis()->GetBinCenter(i);
		
		double tmp_E_beam_targ=targ.beam_e_centre(beam_Z,beam_A,E_beam,TVector3(0.0,0.0,1.0),fraction);
		double tmp_P_beam_targ=get_rel_mom(tmp_E_beam_targ,beam_mass);
		double tmp_beta_CoM=get_com_beta(tmp_P_beam_targ,beam_mass,targ_mass);
		double tmp_KE_0_tot_CoM=get_com_KE(tmp_P_beam_targ,beam_mass,targ_mass);	
		double tmp_P_0_CoM=momentum_beam_P_CoM(tmp_P_beam_targ,beam_mass,targ_mass);
		
		dist_target_KE_beam.SetPoint(dist_target_KE_beam.GetN(),fraction,tmp_E_beam_targ);	
		dist_target_KE_0.SetPoint(dist_target_KE_0.GetN(),fraction,tmp_KE_0_tot_CoM);		
		dist_target_beta.SetPoint(dist_target_beta.GetN(),fraction,tmp_beta_CoM);		
		dist_target_P_0.SetPoint(dist_target_P_0.GetN(),fraction,tmp_P_0_CoM);	
		
		double tmp_reco_E_star=reco_E_star; 
		if(targ_fusion){
			tmp_reco_E_star=tmp_KE_0_tot_CoM+QVal;
			if(tmp_reco_E_star<0)tmp_reco_E_star=0;	
			dist_target_E_star.SetPoint(dist_target_E_star.GetN(),fraction,tmp_reco_E_star);
		}
		
		double tmp_KE_1_tot_CoM=tmp_KE_0_tot_CoM+QVal-tmp_reco_E_star;		
		if(tmp_KE_1_tot_CoM<0)tmp_KE_1_tot_CoM=0;
		dist_target_KE_1.SetPoint(dist_target_KE_1.GetN(),fraction,tmp_KE_1_tot_CoM);

		double tmp_P_1_CoM=momentum_energysplit_CoM(tmp_KE_1_tot_CoM,(reco_mass+(tmp_reco_E_star/jam_phys_amu)),ejec_mass);		
		dist_target_P_1.SetPoint(dist_target_P_1.GetN(),fraction,tmp_P_1_CoM);		
		
		double XXX;
		switch ( dist_target_var ){
			case 1:	XXX=tmp_P_0_CoM;break;
			case 2:	XXX=tmp_beta_CoM;break;
			case 3:	XXX=tmp_reco_E_star;break;
			default:XXX=tmp_KE_0_tot_CoM;
		}
		dist_target_hist.SetBinContent(i,dist_target.Eval(XXX));		
	}
	
	//dist_target_hist.Scale(1/dist_target_hist.Integral());
	dist_target_hist.ComputeIntegral();	
	
	thick_target_fraction_0=dist_target_hist.GetMean();

	thick_beta_CoM=dist_target_beta.Eval(thick_target_fraction_0);
	thick_KE_0_tot_CoM=dist_target_KE_0.Eval(thick_target_fraction_0);
	thick_P_0_CoM=dist_target_P_0.Eval(thick_target_fraction_0);
	if(targ_fusion){thick_reco_E_star=dist_target_E_star.Eval(thick_target_fraction_0);}	
	thick_KE_1_tot_CoM=dist_target_KE_1.Eval(thick_target_fraction_0);
	thick_P_1_CoM=dist_target_P_1.Eval(thick_target_fraction_0);
	this->null_target_pos();
}

void exp_core::gen_random_target_pos(){
	current_target_fraction=dist_target_hist.GetRandom();
	decay_target_fraction=current_target_fraction;

	beta_CoM=dist_target_beta.Eval(current_target_fraction);
	KE_0_tot_CoM=dist_target_KE_0.Eval(current_target_fraction);
	P_0_CoM=dist_target_P_0.Eval(current_target_fraction);
	if(targ_fusion){
		reco_E_star=dist_target_E_star.Eval(current_target_fraction);
	}	
	KE_1_tot_CoM=dist_target_KE_1.Eval(current_target_fraction);
	P_1_CoM=dist_target_P_1.Eval(current_target_fraction);
}

void exp_core::null_target_pos(){if(thick_target_calcs){
	current_target_fraction=thick_target_fraction_0;beta_CoM=thick_beta_CoM;KE_0_tot_CoM=thick_KE_0_tot_CoM;P_0_CoM=thick_P_0_CoM;
	if(targ_fusion){reco_E_star=thick_reco_E_star;}
	KE_1_tot_CoM=thick_KE_1_tot_CoM;P_1_CoM=thick_P_1_CoM;
}}	

void exp_core::lab_loretz_targ_exit(int part,double target_fraction){
	if(part>=0&&part<4){
		if(*part_M[part]<=0)return;//gamma escape
		
		TVector3 Pvec(0,0,0);
		double KER=get_KE(lor_point[part]);
		if(KER>0.001){	
			Pvec=lor_point[part]->Vect();
			KER=targ.particle_e_exit(*part_Z[part],*part_A[part],KER,Pvec,target_fraction);	
			if(KER>0.001){
				Pvec.SetMag(get_rel_mom(KER,*part_M[part]));
			}else{Pvec.SetXYZ(0,0,0);}
		}	
		*lor_point[part]=make_lorentzvec(Pvec,*part_M[part]);
	}
}

