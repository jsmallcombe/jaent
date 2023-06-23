
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
void exp_core::basic_make_single(bool obstructions,int multiR,int multiE,int multiA,int multiB){
	// C++11CODE //vector< int > parhit={0,0,0,0};
	// C++11CODE //vector< int > parmin={multiR,multiE,multiA,multiB};
	vector< int > parhit(4,0);
	vector< int > parmin(4,0);parmin[0]=multiR;parmin[1]=multiE;parmin[2]=multiA;parmin[3]=multiB;

	int sum=0;
	bool matched=false;
	while(!matched&&sum<10000){
		int valid=0;
		this->gen_event();
		this->det_check_hits_all(obstructions);
		
		for(int j=0;j<4;j++){
			int M=det_hits[j].size();
			if(M>=parmin[j]){
				if(M>0)parhit[j]++;
				valid++;
			}
		}
				
		if(valid==4){
			this->addhitsdet();
			matched=true;
		}
		sum++;
	}
}


void exp_core::basic_hit_count(int reps,bool obstructions,int multiR,int multiE,int multiA,int multiB){
	// C++11CODE //vector< int > parhit={0,0,0,0};
	// C++11CODE //vector< int > parmin={multiR,multiE,multiA,multiB};
	vector< int > parhit(4,0);
	vector< int > parmin(4,0);parmin[0]=multiR;parmin[1]=multiE;parmin[2]=multiA;parmin[3]=multiB;

	int validsum=0;
	//decay_events=true;
	for(int i=1;i<=reps;i++){
		int mult=0;
		int valid=0;
		this->gen_event();
		this->det_check_hits_all(obstructions);
		
		for(int j=0;j<4;j++){
			int M=det_hits[j].size();
			if(M>0)mult++;
			if(M>=parmin[j]){
				if(M>0)parhit[j]++;
				valid++;
			}
		}
				
		if(valid==4&&mult>0){this->addhitsdet();validsum++;}
		
		if(i==reps||i==100000||i==1000000||i==10000000||i==100000000||i==1000000000){
			cout<<endl<<endl<<"Total Reaction Events "<<i;
			for(int j=0;j<4;j++){
				if(parhit[j]>0)cout<<endl<<part_names[j]<<" "<<parhit[j]<<" ("<<(double)parhit[j]/(double)i<<") "<<flush;
			}
			cout<<endl<<"Full Events "<<validsum<<" ("<<(double)validsum/(double)i<<") "<<flush;
		}
	}
}

double exp_core::basic_det_hit_frac(int reps,int det,bool obstructions){
	int A=0;
	for(int i=1;i<=reps;i++){
		this->gen_event();
		this->det_check_hits_all(obstructions);
		
		for(int j=0;j<4;j++)
			if(valid_dets[j][det])
				for(int k=0;(unsigned)k<det_hits[j].size();k++)if(det_hits[j][k]==det){A++;j=4;break;}
	}
	return (double)A/(double)reps;
}

double exp_core::basic_det_hit_multi(int reps,int det,bool obstructions){
	int A=0;
	int M=0;
	for(int i=1;i<=reps;i++){
		int B=0;
		this->gen_event();
		this->det_check_hits_all(obstructions);
		
		for(int j=0;j<4;j++)
			if(valid_dets[j][det])
				for(int k=0;(unsigned)k<det_hits[j].size();k++)if(det_hits[j][k]==det){M++;B++;}
		if(B>0)A++;
	}
	return (double)M/(double)reps;
}

void exp_core::basic_decay_check(int reps){
	int s=implant_mode;
	implant_mode=0;
	int a=0,b=0,c=0;
	for(int i=1;i<=reps;i++){
		this->gen_event();
		this->det_check_hits_all(true);
		
		double decay_length =(offset_decay-offset_primary).Mag();
		if(decay_length<0.001)a++;
		else{
			if(stop_distance[0]<=decay_length&&stop_distance[0]>=0)b++;
			else c++;
		}
			
		if(i==reps||i==100000||i==1000000||i==10000000||i==100000000||i==1000000000){
			cout<<endl<<endl<<"Total Reaction Events "<<i;
			cout<<endl<<"Decays in target "<<a<<" ("<<(double)a/(double)i<<") "<<flush;
			cout<<endl<<"Decays after implant "<<b<<" ("<<(double)b/(double)i<<") "<<flush;
			cout<<endl<<"Decays in vacuum "<<c<<" ("<<(double)c/(double)i<<") "<<flush;
		}	
	}
	implant_mode=s;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

//Master function

void exp_core::gen_event(){	
	if(use_offsets&&decay_events){offset_decay.SetXYZ(0,0,0);if(stopper_seperation>0)offset_stopper.SetXYZ(0,0,0);}
// 	offset_primary.SetXYZ(0,0,0);offset_decay
	
	//IF THICK TARGET GENERATE A HIT POSITION IN TARGET AND LOAD ENERGIES
	if(thick_target_calcs)if(targ.GetThickness()>0)this->gen_random_target_pos();

	//Generate the angles and lorentz vectors of the primary recoil+ejectile
	this->gen_primary();

	//Generate decay of the recoil
	if(decay_events)this->gen_decay();
		
	// Now have 2-3 Lorentz vectors boosted to lab frame (but maybe still in target)
	//
	//  Can now get all lab information from TLorentze vectors + offset_decay + offset_primary
	// 		lorH.Phi();
	//		lorH.Theta();
	//		lorH.Vect();
	//		P=lorH.Vect().Mag;
	//		mass=lorH.M()/jam_phys_amu;
	//	BUT IT WILL BE THE NASTY ROOT PHI SO DONT FORGET THAT
	//
	
	//If target energy losses have been requested
	//calc energy loss going out of the target
	if(thick_target_calcs)if(targ.GetThickness()>0)this->gen_Eloss();
		
// 	if(targ.GetThickness()>0&&thick_target_calcs)this->null_target_pos();//reset beam energies etc to nominal values out of paranoia

	if(use_offsets&&decay_events)offset_decay+=offset_primary;//important else decays relative to zero!
	//we only do it at the end of event generation though because relative positions needed during
}

///////////////////////////

void exp_core::gen_primary(){

	TVector3 tmp_vec=TVector3(0.0,0.0,0.0);
	if(!targ_fusion&&!fusion_evaporation_simplify){//save on pointless calculation

		//Generate the in frame back-to-back tragjectories
		//CAN HAVE NON UNIFORM
		//Functions set globals thetaM,phiM which are ejectile (beam like) angles.
		switch ( simm_dist ){
			case 1: this->gen_ruth(thetaM,phiM);break;
			case 2: this->gen_zeroid(thetaM,phiM);break;
			case 3: this->gen_manual_primary(thetaM,phiM);break;
			default: this->gen_uniform(thetaM,phiM);
		}
		
		if(reverse_primary)thetaM=pi-thetaM;
		
		tmp_vec.SetMagThetaPhi(P_1_CoM,thetaM,phiM);
		
		lorE=make_lorentzvec(tmp_vec,ejec_mass);
		lorE.Boost(0,0,beta_CoM);		
	}
	lorR=make_lorentzvec(-tmp_vec,(reco_mass+(reco_E_star/jam_phys_amu)));
	lorR.Boost(0,0,beta_CoM);
	
}

/////////////////////////////////

void exp_core::gen_decay(){
	
	TVector3 labboost;
	
	//The velocity of the recoil in the lab frame directly following primary reaction
	if(targ_fusion){
		labboost=TVector3(0,0,beta_CoM);
		post_beta=beta_CoM;
	}else{
		labboost=lorR.BoostVector();
		post_beta=labboost.Mag();
	}
	
	//If the recoil moves before it decays
	//Calculate where it is and how much post_beta decays
	if(lifetime>0&&post_beta>0){
		//generate a lifetime
		double time_ns_CoM=rand.Exp(lifetime);
		
		//calculate target slowing/stopping
		if(thick_target_calcs&&targ.GetThickness()>0){
			decay_target_fraction=current_target_fraction; //this line might not be needed
			//after this call decay_target_
			post_beta=this->lifetime_track(0,time_ns_CoM,decay_target_fraction);
			labboost.SetMag(post_beta);
		}
		//Now both time and boost have been adjusted for target. If needed.
		
		//if still going after target, offset
		if(time_ns_CoM>0&&post_beta>0){
			double remaining_lab_ns=time_ns_CoM*get_gamma(post_beta);
			offset_decay=labboost;//gets direction
			offset_decay.SetMag(jam_phys_speed_c_mm_ns*post_beta*remaining_lab_ns);	
		}
		
		if(stopper_seperation>0){offset_decay+=offset_stopper;}		
	}
	// need to consider hitting things

	
	//Generate the selected kind of decay
	switch ( decay_type ){
		case 1://currently nothing to do for a gamma ray
			break;
		case 2:
			//sets beta AND recoil to have same momentum, acceptable simplifications
			P_decay_recoil_CoM=beta_dist_hist.GetRandom();
			break;
		case 3:	//currently nothing to do for a ICE	
			break;
		case 4:	//currently nothing to do for a alpha decay	
			break;	
		default://fission
			this->gen_fission();
			P_decay_recoil_CoM=momentum_energysplit_CoM(KE_decay_tot_CoM,mass_decay1,mass_decay2);
	}

	//Generate the back-to-back angles in the decay frame
	//Angles are relative to beam currently, but lorR is avaible for rotation
	//CAN ADD NON UNIFORM HERE
	//this->gen_uniform();
	phiS=pi*rand.Uniform(0,2);
	thetaS=acos(rand.Uniform(-1,1));
	
	TVector3 tmp_vec;		
	tmp_vec.SetMagThetaPhi(P_decay_recoil_CoM,thetaS,phiS);
	lorA=make_lorentzvec(tmp_vec,mass_decay1+(decay_recoil_E_star/jam_phys_amu));
	lorB=make_lorentzvec(-tmp_vec,mass_decay2);

// 	if(targ_fusion){	
// 		//boost back to the lab frame
// 		lorA.Boost(0,0,post_beta);
// 		lorB.Boost(0,0,beta_CoM);				
// 	}
	if(post_beta>0){
		lorA.Boost(labboost);
		lorB.Boost(labboost);
	}	
}

///////////////////////////////


void exp_core::gen_Eloss(){
	
	//recoil and ejectile attempt to leave target
	this->lab_loretz_targ_exit(0,current_target_fraction);
	if(!targ_fusion){	
		this->lab_loretz_targ_exit(1,current_target_fraction);
	}
	
	
	if(decay_events){
		double m=0,a=0;
			//determine the distance from the surface of the target
		if(use_offsets&&(decay_target_fraction>=2||decay_target_fraction<=-1)){
			m=offset_decay.Mag();
			a=targ.GetNormal().Angle(offset_decay);
			m*=abs(cos(a));
		}
		
		for(int i=2;i<4;i++){//loop for decay daughter and particle
			double d=0;
			
			//if it's left the target calculate how close it gets to the origin
			if(m>0){
				TVector3 lv=lor_point[i]->Vect();
				double n=lv.Mag();
				double b=targ.GetNormal().Angle(lv);					
				
				if(n>0&&((a<pi*0.5)^(b<pi*0.5))){
					lv.SetMag(m/abs(cos(b)));
					d=(offset_decay+lv).Mag();
					
				}else{d=1000;}
			}
		
			//if it never left OR passes within 10mm then it's going through the target
			if(d<10.0)this->lab_loretz_targ_exit(i,decay_target_fraction);
			
			if(decay_type==1||decay_type==2||decay_type==3) i++;//no energy loss calcs for gammas nor implemented for beta/ICE
		}
	}
}
