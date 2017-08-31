
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

void exp_core::set_target_primary_offset(TVector3 offin){
	offset_primary=offin;
	if(offin.Mag()>0)use_offsets=true;
	else if(lifetime<=0) use_offsets=false;
}
void exp_core::set_target_primary_offset(double X,double Y,double Z){
	this->set_target_primary_offset(TVector3(X,Y,Z));
}

void exp_core::set_lifetime_ns(double lifetime_in,double a,double b){
	if(lifetime_in>0){
		lifetime=abs(lifetime_in);
		use_offsets=true;
		if(a>0)density_targ=a;
		if(b>0)density_back=b;
		this->print_target();
	}else{
		lifetime=0;
		if(offset_primary.Mag()<=0)use_offsets=false;
	}
}

void exp_core::set_halflife_ns(double halflife_in,double a,double b){this->set_lifetime_ns(halflife_in/log(2),a,b);};
		
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
	

double exp_core::get_range(int i,bool back,double mev){//output in mg
	if(mev>0){
		if(back)return targ.GetRangeBack(*part_Z[i],*part_A[i],mev);
		else return targ.GetRange(*part_Z[i],*part_A[i],mev);
	}
	return 0.0;
}

double exp_core::passage(int i,bool back,double mev,double mgcm2){
	if(mev>0){
		if(back)return targ.passB(*part_Z[i],*part_A[i],mgcm2,mev);
		else return targ.pass(*part_Z[i],*part_A[i],mgcm2,mev);
	}
	return 0.0;
}


 
double exp_core::lifetime_track(int i,double& ns,double& fraction_in){
	//Given the particle ID, the fraction into the target and time on the clock
	//ID gives access to initial direction, mass and momentum.
	//This function then calculates target passage
	//Returns velocity either and the moment the clock runs out OR at point of leaving target
	//Updates fraction_in	
	
	//calculate initial parameters
	double KE_0=get_KE(lor_point[i]);	//KE in lab frame at point of interaction (fraction_in==current_target_fraction for this event initially)
	double beta_0=lor_point[i]->Beta();	//beta of parent nucleus
	double mass= *part_M[i];		//mass of parent nucleus in amu
	TVector3 traj=lor_point[i]->Vect();	//momentum vector lab frame
	
	double sep=0;//place holder for a stopper separation, in mm whereas targets in um

	//effective thickness factor from angle > 1
	double eff_thick=target_effective(targ.GetNormal(),traj,1.0);
	
	//some logic gate about the relative layers and trajectories
	bool beamward=false,backingward=false,b_start=false;
	if(abs(targ.GetFormal().Angle(traj))<pi/2)beamward=true;//going out the "front" i.e. frac goes up
	if(targ.GetBackThickness()>0 && (targ.GetNormal().Angle(traj)>pi/2))backingward=true;//There IS a backing and trajectory "towards" backing 
	if(abs(fraction_in-0.5)>0.5&&targ.GetBackThickness()>0)b_start=true;//If (for some unknown reason) starting in the backing
	bool two_layer=false;
	if(b_start^backingward)two_layer=true;// If two layers will be traversed in this trajectory and starting position

	if(stopper_seperation>0){ //If there is a stopper gap (stopper concept quite limited)
		sep=stopper_seperation*eff_thick;// Length of path between stoppers
		
		// Set the offset (relative to target position) of the particle if it exits through the stopper
		if(backingward)offset_stopper=traj.Unit()*sep;
		else offset_stopper=traj.Unit()*-sep;
	}
	
	//fraction_in defined as such
	//-1 -> 0 upstream backing/stopper
	// 0 -> 1 target
	// 1 -> 2 downstream backing/stopper
	//In the following we convert fraction_in to a fraction of the current layer (frac1 0->1) between the current position and exit of the layer (up or downstream, depending on input trajectory)
	
	// Data on the thickness in um and mg/cm2 of the layers being traversed
	// Partial thickness for the one we are in and total for the 2nd if it exists
	double frac1=fraction_in;	
	double thick1=targ.GetThickness()*eff_thick;
	double thick2=targ.GetBackThickness()*eff_thick;
	double um1=giveme_um(density_targ,thick1,false);
	double um2=giveme_um(density_back,thick2,false);	
	if(b_start){
		frac1=fraction_in-1.0;//set distance from upstream
		if(fraction_in<0)frac1=fraction_in+1.0;//set distance from upstream

		swap_jd(um1,um2);
		swap_jd(thick1,thick2);
	}
	if(beamward)frac1=1-frac1;//only needed for frac1 as frac2 is always 1!
	um1*=frac1;//only needed for frac1 as frac2 is always 1!

	////////////////////////////////////////// 
	///////// start actual calculations  /////		
	//////////////////////////////////////////	
	
	//calc data for complete target exit

	bool stopped=false;
	double ran1=1,ran2=1;//Range in the material of each layer, expressed as a fraction of the layer (can be larger than 1)
	double stopfrac=0,KE_exit=0,beta_exit=0;

	double KE_betw=this->passage(i,b_start,KE_0,thick1*frac1); //Energy after layer 1
	ran1 = get_range(i,b_start,KE_0)/thick1;	//Range in layer as a fraction of the layer
	if(KE_betw>0&&two_layer){
		KE_exit=this->passage(i,!b_start,KE_betw,thick2);//Energy on exit
		ran2 = get_range(i,!b_start,KE_betw)/thick2;//Range in layer as a fraction of the layer
	}else{KE_exit=KE_betw;}//Energy on exit
	if(KE_exit>0) beta_exit=get_beta_KE(KE_exit,mass);//If it exits, what beta

	//The above is only concerned with energy so a vacuum gap between the two in the case of a stopper is of negligible difference

	double x_0=0;//Total distance to stop (or exit) along current trajectory, in um
	
	if(ran1<frac1){//Stopped in layer 1
		stopped=true;
		stopfrac=ran1;
		um1*=(ran1/frac1);//shrink the range down
		um2=0;
		KE_betw=0.0;
		two_layer=false;
		KE_exit=0.0;
		beta_exit=0.0;
		ran2=0.0;
	}else{ran1=frac1;}
	x_0+=um1;
	

	if(two_layer){//Stopped in layer 2
		if(ran2<1){
			stopped=true;
			stopfrac=ran2+frac1;
			um2*=ran2;
			KE_exit=0.0;
			beta_exit=0.0;
		}else{ran2=1;}
		x_0+=um2;
	}
	//ran1 and ran2 are not = frac1 and 1 or the stopping frac if stopped
	//um1 and um2 are updated if needed
	

	//////////////////////// 
	///////// start slowing stuff
	////////////////////////
	
	//Very rough calc how long it takes to stop (or exit)
	//This assumes same deceleration in both, which it shouldn't
	double rough_time_ns=(x_0/((beta_0+beta_exit)*0.5*jam_phys_speed_c_mm_ns*1000));//lab time (trapezoid)
	rough_time_ns/=0.5*(get_gamma(beta_0)+get_gamma(beta_exit));//frame time (dilate average)
	
	// We calculated the time to traverse the stopper BEFORE the time to traverse the gap and reach it
	// But as we are only interested in the total, the order we sum them is inconsequential 
	if(two_layer&&sep>0){ //if there is a stopper gap and it is traversed  (stopper concept quite limited)
		double beta_between=get_beta_KE(KE_betw,mass);
		rough_time_ns+=(sep/(jam_phys_speed_c_mm_ns*beta_between))/get_gamma(beta_between);		
	}
	
	//if exit with time to spare	
	if(beta_exit>0&&!stopped){
		if(rough_time_ns*5<ns){//target interaction finished well before time out
			ns=ns-rough_time_ns;	

			if(backingward)fraction_in=-1+(3*beamward);
			else fraction_in=1.0*beamward;		

			if(stopper_seperation>0&&!backingward)offset_stopper*=0;//didnt fly out stopper	
			
			return beta_exit;
		}
	}
	
	//if stopped with time to spare	
	if(stopped){
		if(rough_time_ns*2<ns){//target interaction finished well before time out
			ns=0.0;	

			if(beamward) fraction_in+=stopfrac;
			else fraction_in-=stopfrac;	
			
			if(fraction_in>=0&&fraction_in<=1)offset_stopper*=0;//stopped in target	

			return 0.0;
		}
	}	
	

// 	//IF we have reached this point, EITHER stops shortly before decay OR decays in/near target
// 	//Do two careful loop to see if realtime expires before exit OR stop
// 	//The loops run until particle exit if it will exit OR particle stop if it will stop
// 	
	//starting conditions into loop variables
	double KE_c=KE_0;
	double beta_c=beta_0;
	double t_t=0.0;

// 	//do 10 steps through the fraction of the first layer we must traverse
	for(int j=0;j<10;j++){
		double KE_l=this->passage(i,b_start,KE_c,thick1*0.1*ran1);
		
		//time taken to cross this section		
		double beta_l=get_beta_KE(KE_l,mass);
		double t_l=((um1*0.1)/((beta_c+beta_l)*0.5*jam_phys_speed_c_mm_ns*1000));
		t_l/=0.5*(get_gamma(beta_c)+get_gamma(beta_l));	
		
		// if decayed in that 10th
		if(t_l+t_t>=ns){
			double partway=abs(ns-t_t)/t_l;
			double KEF=KE_c+((KE_c-KE_l)*partway);
			partway=0.1*((double)j+partway)*ran1;

			if(beamward) fraction_in+=partway;
			else fraction_in-=partway;

			ns=0.0;//remaining after target
			if(fraction_in>=0&&fraction_in<=1)offset_stopper*=0;//stopped in target	
			return get_beta_KE(KEF,mass);//beta
		}			
			
		//prepare for next loop
		t_t+=t_l;
		KE_c=KE_l;
		beta_c=beta_l;
		if(KE_c==0){j=10;break;}
	}
	
// 	//gap between layers "stopper"
	if(two_layer&&sep>0){ //if there is a stopper gap and it is traversed  (stopper concept quite limited)
		if(KE_c>0){
			double t_g=(sep/(jam_phys_speed_c_mm_ns*beta_c))/get_gamma(beta_c); //does it decay in flight?
			if(t_g+t_t>=ns){//decayed while passing gap

				double ftr=(ns-t_t)/t_g; //shrink the offset vector
				if(backingward) offset_stopper*=ftr;
				else offset_stopper*=(1-ftr);
				
				if(beamward) fraction_in=1;
				else fraction_in=0;

				ns=0.0;//remaining after target			
				return beta_c;//beta
			}
			t_t+=t_g;//no decay
		}
		
		if(!backingward)offset_stopper*=0;//didnt fly out stopper	
	}

// 	//do 10 steps through the fraction of the second layer we must traverse	
	// ran2 should = 1 if particle doesnt stop
	if(two_layer){
		KE_c=KE_betw;

		//do 10 steps through the fraction of the second layer we must traverse
		for(int j=0;j<10;j++){
			double KE_l=this->passage(i,!b_start,KE_c,thick2*0.1*ran2);

			//time taken to cross this section		
			double beta_l=get_beta_KE(KE_l,mass);
			double t_l=((um2*0.1)/((beta_c+beta_l)*0.5*jam_phys_speed_c_mm_ns*1000));
			t_l/=0.5*(get_gamma(beta_c)+get_gamma(beta_l));

			//if decayed in that 10th
			if(t_l+t_t>=ns){
				double partway=abs(ns-t_t)/t_l;
				double KEF=KE_c+((KE_c-KE_l)*partway);
				partway=0.1*((double)j+partway)*ran2;

				if(beamward) fraction_in+=(partway+ran1);
				else fraction_in-=(partway+ran1);

				ns=0.0;//remaining after target			
				return get_beta_KE(KEF,mass);//beta
			}			
				
			//prepare for next loop
			t_t+=t_l;
			KE_c=KE_l;
			beta_c=beta_l;
			if(KE_c==0){j=10;break;}
		}
	}
	
	
	// IF it stops in the target and doesnt decay before this point
	if(stopped){
		ns=0.0;	

		if(beamward) fraction_in+=stopfrac;
		else fraction_in-=stopfrac;	

		return 0.0;
	}	
	
	
	//if it reached this point it has exited without decaying	
	if(backingward)fraction_in=-1+(3*beamward);
	else fraction_in=1.0*beamward;	
	ns-=t_t;//remaining after target			
	return beta_exit;//beta	
}
 
