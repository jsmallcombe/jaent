
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

void exp_core::set_fission(double mass_ratio_in,double delta_m_in,double delta_TKE_in,int TKECentral_in){
	decay_events=true;
	decay_type=0;
	decay_recoil_E_star=0;
	ratioAZ = (double)reco_Z/(double)reco_A;
	if(TKECentral_in==0)TKECentral=true;else TKECentral=false;
	
	mass_ratio=mass_ratio_in;
	if(mass_ratio>2)mass_ratio=1;
	if(mass_ratio>1)mass_ratio=2-mass_ratio;
	delta_m=delta_m_in;
	delta_TKE=delta_TKE_in;
	
	mass_decay1=(reco_mass*0.5*mass_ratio);
	mass_decay2=reco_mass-mass_decay1;	
	
	decay_A_A=round(mass_decay1);
	decay_A_B=reco_A-decay_A_A;	
	decay_Z_A=(mass_decay1*ratioAZ)+0.5;//lazy int rounding
	decay_Z_B=reco_Z-decay_Z_A;	
	
	if(TKECentral_in>0)
		frag_bar_max=TKECentral_in;
		else
		frag_bar_max=viola_TKE(reco_A,reco_Z);
	frag_bar_max/=classical_barrier(reco_A*0.5,reco_Z*0.5,reco_A*0.5,reco_Z*0.5);
	TKE_fiss_zero=frag_bar_max*classical_barrier(decay_A_B,decay_Z_A,decay_A_B,decay_Z_B);
	KE_decay_tot_CoM=TKE_fiss_zero;

	P_decay_recoil_CoM=momentum_energysplit_CoM(KE_decay_tot_CoM,mass_decay1,mass_decay2);
}


TH2D* exp_core::fission_fragment_dist(){
	TH2D* ret = new TH2D("fragment_dist","fragment_dist",int(TKE_fiss_zero*delta_TKE*12),TKE_fiss_zero-TKE_fiss_zero*delta_TKE*6,TKE_fiss_zero+TKE_fiss_zero*delta_TKE*6,int((reco_mass*0.5*(1-mass_ratio)+delta_m*6)*2),reco_mass*0.5*mass_ratio-delta_m*6,reco_mass*0.5*(2-mass_ratio)+delta_m*6);
	ret->GetXaxis()->SetTitle("TKE (MeV)");
	ret->GetYaxis()->SetTitle("Fragment Mass (amu)");
	
	for(int i=0;i<500000;i++){
		this->gen_fission();
		ret->Fill(KE_decay_tot_CoM,mass_decay1);
		ret->Fill(KE_decay_tot_CoM,mass_decay2);	
	}
	return ret;
}
TH2D* exp_core::fission_fragment_dist_detect(bool obstructions){
	TH2D* ret = new TH2D("fragment_dist_detect","fragment_dist_detect",int(TKE_fiss_zero*delta_TKE*12),TKE_fiss_zero-TKE_fiss_zero*delta_TKE*6,TKE_fiss_zero+TKE_fiss_zero*delta_TKE*6,int((reco_mass*0.5*(1-mass_ratio)+delta_m*6)*2),reco_mass*0.5*mass_ratio-delta_m*6,reco_mass*0.5*(2-mass_ratio)+delta_m*6);
	ret->GetXaxis()->SetTitle("TKE (MeV)");
	ret->GetYaxis()->SetTitle("Fragment Mass (amu)");
	
	int points=0;
	int reps=0;
	while(points<500000){
		reps++;
		this->gen_event();
		this->det_check_hits_all(obstructions);	
		
		if(det_hits[2].size()>0&&det_hits[3].size()>0){
			points++;
			ret->Fill(KE_decay_tot_CoM,mass_decay1);
			ret->Fill(KE_decay_tot_CoM,mass_decay2);
		}
		
		if(reps>1000000&&((double)points/(double)reps<0.001))break;
	}
	return ret;
}

TH2D* exp_core::fission_fragment_angles(int N){
	TH2D* ret = new TH2D("fragment_angles","fragment_angles",200,1,4,100,-0.5,0.5);
	ret->GetXaxis()->SetTitle("dTheta");
	ret->GetYaxis()->SetTitle("dPhi");
	
	for(int i=0;i<N;i++){
		this->gen_event();
		
		double a=lorA.Vect().Theta();
		double b=lorB.Vect().Theta();
		
		if(a>0&&b>0){
			double d_theta=a+b;
			double d_phi=angledifference_signed(lorA.Vect().Phi(),lorB.Vect().Phi()-pi);
			
			ret->Fill(rand.Gaus(d_theta,0.02),rand.Gaus(d_phi,0.02));
		}
	}
	return ret;
}
TH2D* exp_core::fission_fragment_angles_detect(bool obstructions,int N){
	TH2D* ret = new TH2D("fragment_angles_detect","fragment_angles_detect",200,1,4,100,-0.5,0.5);
	ret->GetXaxis()->SetTitle("dTheta");
	ret->GetYaxis()->SetTitle("dPhi");
	
	int points=0;
	int reps=0;
	while(points<N){
		reps++;
		this->gen_event();
		this->det_check_hits_all(obstructions);	
		
		if(det_hits[2].size()>0&&det_hits[3].size()>0){
			points++;
		
			double d_theta,d_phi;
			if(use_offsets&&offset_decay.Mag()>0){
				TVector3 vectA=this->path_vec_calc(2,det_hits[2][0]);	
				TVector3 vectB=this->path_vec_calc(3,det_hits[3][0]);	
				vectA+=offset_decay;	
				vectB+=offset_decay;
				
				d_theta=vectA.Theta()+vectB.Theta();
				d_phi=angledifference_signed(vectA.Phi(),vectB.Phi()-pi);
			}else{
				d_theta=lorA.Vect().Theta()+lorB.Vect().Theta();
				d_phi=angledifference_signed(lorA.Vect().Phi(),lorB.Vect().Phi()-pi);
			}
			
			ret->Fill(rand.Gaus(d_theta,0.02),rand.Gaus(d_phi,0.02));	
		}
		if(reps>1000000&&((double)points/(double)reps<0.001))break;
	}
	return ret;
}
TH3D* exp_core::fission_fragment_transfer_angles(bool obstructions){
	TH3D* ret = new TH3D("fission_fragment_transfer_angles","fission_fragment_transfer_angles",200,1,4,100,-0.5,0.5,100,0,pi);
	ret->GetXaxis()->SetTitle("dTheta");
	ret->GetYaxis()->SetTitle("dPhi");
	
	int points=0;
	int reps=0;
	while(points<500000){
		reps++;
		this->gen_event();
		this->det_check_hits_all(obstructions);	
		
		if(det_hits[2].size()>0&&det_hits[3].size()>0){
			points++;
		
			double d_theta=lorA.Vect().Theta()+lorB.Vect().Theta();
			double d_phi=angledifference_signed(lorA.Vect().Phi(),lorB.Vect().Phi()-pi);
			d_phi=rand.Gaus(d_phi,0.01);	
			
			ret->Fill(d_theta,d_phi,lorE.Vect().Theta());	
		}
		if(reps>1000000&&((double)points/(double)reps<0.001))break;
	}
	return ret;
}


TH1D* exp_core::fission_fragment_dT(int detA,int detB,bool obstructions){
	TH1D* ret = new TH1D("fragment_dT","fragment_dT",400,-100,100);
	ret->GetXaxis()->SetTitle("ns");
	ret->GetYaxis()->SetTitle("counts");
	
	int points=0;
	int reps=0;
	while(points<100000){
		reps++;
		this->gen_event();
		this->det_check_hits_all(obstructions);	
		
		if(det_hits[2].size()>0&&det_hits[3].size()>0){
			if((det_hits[2][0]==detA||det_hits[2][0]==detB)&&(det_hits[3][0]==detA||det_hits[3][0]==detB)){//if selected detectors are hit
				points++;
				
				double t1=(detectors[det_hits[2][0]].hit_pos_vec(lor_point[2],offsets[2]).Mag())/(lor_point[2]->Beta()*jam_phys_speed_c*0.000001);
				double t2=(detectors[det_hits[3][0]].hit_pos_vec(lor_point[3],offsets[3]).Mag())/(lor_point[3]->Beta()*jam_phys_speed_c*0.000001);
				if(det_hits[2][0]==detA)ret->Fill(t1-t2);
				else ret->Fill(t2-t1);
			}
		}//cout<<reps<<" "<<points<<endl;
		if(reps>1000000&&((double)points/(double)reps<0.001))break;
	}
	return ret;
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

// void exp_core::(){
// }
void exp_core::gen_fission(){
	KE_decay_tot_CoM=TKE_fiss_zero;
	if(delta_m>0){
		mass_decay1=rand.Gaus(reco_mass*0.5*mass_ratio,delta_m);
		if(rand.Integer(2)==0)mass_decay1=reco_mass-mass_decay1;//randomise which is lightest
		mass_decay2=reco_mass-mass_decay1;	
		decay_A_A=round(mass_decay1);
		decay_A_B=reco_A-decay_A_A;	
		decay_Z_A=(mass_decay1*ratioAZ)+0.5;//lazy int rounding
		decay_Z_B=reco_Z-decay_Z_A;
		if(!TKECentral)KE_decay_tot_CoM=frag_bar_max*classical_barrier(decay_A_A,decay_Z_A,decay_A_B,decay_Z_B);
	}	
	if(delta_TKE>0){
		KE_decay_tot_CoM=rand.Gaus(KE_decay_tot_CoM,KE_decay_tot_CoM*delta_TKE);

	}
	if(delta_TKE>0||delta_m>0){
		P_decay_recoil_CoM=momentum_energysplit_CoM(KE_decay_tot_CoM,mass_decay1,mass_decay2);
	}
}