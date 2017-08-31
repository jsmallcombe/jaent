
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

void exp_core::print_reaction(){
	if(targ.GetThickness()>0&&thick_target_calcs)this->null_target_pos();
	
	if(targ_fusion){
		cout<<endl<<endl<<"  "<<this->channelnames(beam_Z,beam_A);
		cout<<" + "<<this->channelnames(targ_Z,targ_A);
		cout<<" -> "<<this->channelnames(reco_Z,reco_A);
	}else{
		cout<<endl<<endl<<"  "<<this->channelnames(targ_Z,targ_A);
		cout<<"("<<this->channelnames(beam_Z,beam_A);
		cout<<","<<this->channelnames(ejec_Z,ejec_A);
		cout<<")"<<this->channelnames(reco_Z,reco_A);		
	}
	if(reco_E_star>0)cout<<"*";
	if(decay_events){
		cout<<" -> ";
		if(ejec_A>0)cout<<this->channelnames(ejec_Z,ejec_A)<<" + ";
		cout<<this->channelnames(decay_Z_A,decay_A_A)<<" + ";
		cout<<this->channelnames(decay_Z_B,decay_A_B);
	}
	
	cout<<endl<<"Q = "<<QVal<<" MeV";
	if(!validQ)cout<<" : INVALID Q, INPUT MANUALLY."<<endl;
	if(QMan!=0.0)cout<<" : MANUAL Q.";
	
	if(reco_E_star>0)cout<<"  E* = "<< reco_E_star<<" MeV";
	
	double KEsm=QVal-reco_E_star;
	if(KEsm>0)KEsm=0;
	cout<<"  KE_CoM_min = "<<abs(KEsm)<<" MeV";
	
	cout<<"  Coul_barrier = "<<basic_barrier_in<<" MeV";

	cout<<endl<<endl<<"     E_Beam_min = "<<E_beam_min<<" MeV";
	cout<<endl<<" E_Beam_barrier = "<<E_beam_barrier<<" MeV";
	cout<<endl<<"         E_beam = "<<E_beam<<" MeV";
	cout<<"  Beta = "<<beta_CoM;	
	cout<<"  E_CoM = "<<KE_0_tot_CoM<<" MeV";	
	if(KE_0_tot_CoM+QVal-reco_E_star<0){
		cout<<endl<<"E BEAM BELOW THRESHOLD"; 
	}	
	
	
	if(decay_events){
		if(lifetime>0){
			cout<<endl<<" Lifetime = ";
			if(lifetime>=1E9)cout<<lifetime/1E9<<" s";
			if(lifetime>=1E6&&lifetime<1E9)cout<<lifetime/1E6<<" ms";
			if(lifetime>=1E3&&lifetime<1E6)cout<<lifetime/1E3<<" us";
			if(lifetime>=1&&lifetime<1E3)cout<<lifetime<<" ns";
			if(lifetime>=1E-3&&lifetime<1)cout<<lifetime*1E3<<" ps";
			if(lifetime<1E-3)cout<<lifetime*1E6<<" fs";
		}
		
		switch ( decay_type ){
			case 1:
				if(decay_recoil_E_star<0)cout<<endl<<" INVALID GAMMA ENERGY";
				cout<<endl<<" Ey = "<<P_decay_recoil_CoM*1000<<" keV. E* Final = "<<decay_recoil_E_star;
				break;
			case 2:
			//	cout<<
				break;
			case 3:
			//	cout<<
				break;
			default://fission
				cout<<endl<<" Nominal TKE : "<<TKE_fiss_zero<<" MeV";
		} 
		
	}
	
	cout<<endl;
}


void exp_core::print_decay(){
	cout<<endl<<endl;
	if(decay_type>0&&decay_type<5){
		bool ye=false;//is it just a gamma/electron decay
		if(decay_type==1||decay_type==3)ye=true;
		
		if(reco_E_star>0){
			cout<<setw(6)<<jsigfig(reco_E_star,4)<<"__________"<<endl;
			cout<<setw(17)<<"|"<<endl;
			cout<<setw(18)<<"| ";
			if(ye){
				cout<<jsigfig(reco_E_star-decay_recoil_E_star,4)<<" "<<this->channelnames(decay_Z_B,decay_A_B)<<endl;
				if(decay_recoil_E_star>0){
					cout<<setw(6)<<jsigfig(decay_recoil_E_star,4)<<"__________|"<<endl<<endl<<endl;
				}
			}else{cout<<endl;}
		}
		cout<<setw(6)<<"0"<<"__"<<setw(5)<<this->channelnames(reco_Z,reco_A)<<"___";
		if((reco_E_star>0&&!ye)||(decay_recoil_E_star<=0&&ye))cout<<"|"<<endl;
		else{cout<<endl;}
		
		if(!ye){
			double TQ=(reco_mass-mass_decay1-mass_decay2)*jam_phys_amu;
			double E=TQ+reco_E_star-decay_recoil_E_star;
			
			cout<<setw(17)<<"|"<<setw(18)<<"|"<<endl;
			cout<<setw(17)<<"|"<<setw(6)<<jsigfig(E,4)<<" "<<setw(2)<<this->channelnames(decay_Z_B,decay_A_B)<<setw(9)<<"|"<<setw(5)<<jsigfig(TQ,4)<<endl;
			cout<<setw(17)<<"|"<<setw(18)<<"|"<<endl<<setw(17)<<"|"<<setw(18)<<"|"<<endl<<setw(17)<<"|";
			if(decay_recoil_E_star>0){
				cout<<"__________ "<<setw(5)<<jsigfig(decay_recoil_E_star,4)<<" |"<<endl;
				cout<<setw(35)<<"|"<<endl;
				cout<<setw(35)<<"|"<<endl;
				cout<<setw(17)<<" ";
			}
			cout<<"__"<<setw(5)<<this->channelnames(decay_Z_A,decay_A_A)<<"___ 0     |";
		}
		
	}

	
}






void exp_core::print_detectables(){
	//0 = Recoil, 1 = Ejectile, 2 = decay_A (recoil), 3 = decay_B (eject)
	cout<<endl<<endl<<"               |";
	cout<<setw(5)<<this->channelnames(reco_Z,reco_A)<<"|";
	cout<<setw(5)<<this->channelnames(ejec_Z,ejec_A)<<"|";
	
	if(decay_events&&decay_type==0){
		cout<<"  FF | FF  |";
	}else{
		cout<<setw(5)<<this->channelnames(decay_Z_A,decay_A_A)<<"|";
		cout<<setw(5)<<this->channelnames(decay_Z_B,decay_A_B)<<"|";
	}
	cout<<endl<<"               |Recoi|Eject|Recoi|Decay|";
				
	cout<<endl;
	for(int i=0;i<this->detN();i++){
		cout<<setw(15)<<GetDetName(i)<<"|";
		for(int j=0;j<4;j++){
			if(valid_dets[j][i])cout<<"  *  |";else cout<<"     |";
		}
		cout<<endl;
	}
}


void exp_core::print_doubles()
{
	cout<<endl<<endl<<"Det";
	for(int j=0;(unsigned)j<valid_doubles.size();j++)cout<<setw(2)<<j<<" ";
	cout<<endl;
	for(int i=0;(unsigned)i<valid_doubles.size();i++){
		cout<<setw(2)<<i<<" ";
		for(int j=0;(unsigned)j<valid_doubles[i].size();j++){
			cout<<" "<<valid_doubles[i][j]<<" ";
		}
		cout<<endl;
	}
}

void exp_core::print_target(){
	
	cout<<endl<<" Target "<<this->channelnames(targ_Z,targ_A);
	if(targ.GetThickness()>0){
		cout<<endl<<" Thickness "<<targ.GetThickness()<<" mg/cm2. Density "<<density_targ<<" g/cm3. Thickness "<<giveme_um(density_targ,targ.GetThickness(),false)<<" um"<<flush;
	}
	if(targ.GetBackThickness()>0){
		cout<<endl<<" Backing "<<this->channelnames(targ.bZ(),targ.bA());
		cout<<endl<<" Thickness "<<targ.GetBackThickness()<<" mg/cm2. Density "<<density_back<<" g/cm3. Thickness "<<giveme_um(density_back,targ.GetBackThickness(),false)<<" um"<<flush;
	}cout<<endl;
}


string exp_core::channelnames(int Z,int A){return nuclear_data_ob::channelname(Z,A);}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
