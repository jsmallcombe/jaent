
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



double exp_core::set_uniform(double thetamin,double thetamax){ //input rad,rad
	simm_dist=0;

	if(thetamin>thetamax){double xx=thetamax;thetamax=thetamin;thetamin=xx;}
	if(thetamin<0)thetamin=0;
	if(thetamax>pi)thetamax=pi;

	simm_dist_parA=cos(thetamin);
	simm_dist_parB=cos(thetamax);
	return (cos(thetamin)-cos(thetamax))/2;
}

double exp_core::set_rutherford(double thetamin,double thetamax){ //input rad,rad output in mb
	simm_dist=1;
	this->set_elastic();

	if(thetamin>thetamax){double xx=thetamax;thetamax=thetamin;thetamin=xx;}
	if(thetamin<1E-6)thetamin=1E-6;
	if(thetamax>pi)thetamax=pi;
		
	simm_dist_parA=1/pow(sin(thetamin/2),2);//thetamin
	simm_dist_parB=1/pow(sin(thetamax/2),2);//thetamax

	double cross_easy = rutherford_crosssection(beam_Z,targ_Z,KE_0_tot_CoM,thetamin,thetamax);
	
	//If thick target, cross section is an even sum (simple average) of all the partials summed,
	//Events are then of course distributed through the target weighted by the partial cross-sections
	if(thick_target_calcs)if(targ.GetThickness()>0){
		this->set_target_interaction(3);
		double sumpartcrosssec=0;
		for(int i=0;i<dist_target_KE_0.GetN();i++){
			double x,y;
			dist_target_KE_0.GetPoint(i,x,y);
			sumpartcrosssec+=rutherford_crosssection(beam_Z,targ_Z,y,thetamin,thetamax);
		}sumpartcrosssec/=dist_target_KE_0.GetN();
		return sumpartcrosssec;
	}
	
	return cross_easy;
}

double exp_core::get_rutherford_lab(double thetamin,double thetamax,bool ejec_v,bool reco_v){ //input rad,rad lab output in mb
	this->set_elastic();
	double* a=this->generate_mask(thetamin,thetamax,ejec_v,reco_v);//just to get the cm angles
	double cxr=0;
	for(int i=0;i<4;i++){
		if(a[i]>0&&a[i+4]>0){
			cxr+=this->set_rutherford(a[i],a[i+4]);//we use this because it does thick target stuff
		}
	}
	delete a;
	return cxr;
}

double exp_core::set_rutherford_lab(double thetamin,double thetamax){ //input rad,rad lab output in mb
	double cxr=get_rutherford_lab(thetamin,thetamax);
	this->set_ruthmask(thetamin,thetamax,true,true);
	return cxr;
}

void exp_core::set_zeroid(double sigma){
	simm_dist=2;
	simm_dist_parA=sigma;
}

void exp_core::set_primary_dist(string file){
	ifstream infile;
	infile.open(file);
	if(infile.is_open()){
		double x,y;
		TGraph* bob=new TGraph();
		while(infile>>x>>y){
			bob->SetPoint(bob->GetN(),x,y);
		}
		this->set_primary_dist(bob);
		delete bob;
		infile.close();
		return;
	}
	cout<<endl<<"FILE READ ERROR"<<endl;
}

void exp_core::set_primary_dist(TGraph* data){
	if(data->GetN()>1){
		double a,b,c,d;
		data->Sort();
		data->ComputeRange(a,b,c,d);
		
		double scale=1;
		if((180-c)<abs(c-pi))scale =180/pi;
		
		for(int i=1;i<=manual_primary_dist_hist.GetNbinsX();i++){
			double y=data->Eval(scale*manual_primary_dist_hist.GetXaxis()->GetBinCenter(i));
			manual_primary_dist_hist.SetBinContent(i,y);
		}
		
		manual_primary_dist_hist.Scale(1/manual_primary_dist_hist.Integral());
		manual_primary_dist_hist.ComputeIntegral();
		
		for(int i=1;i<=manual_primary_dist_hist.GetNbinsX();i++){
			manual_primary_dist_store.SetBinContent(i,manual_primary_dist_hist.GetBinContent(i));
		}
		
		simm_dist=3;
		this->reset_mask();
	}
}

void exp_core::set_primary_dist(TF1* data){if(data){
		this->set_primary_dist(data->GetFormula());
}}

void exp_core::set_primary_dist(TFormula* data){if(data){
	TGraph* pass=new TGraph();
	for(int i=1;i<=manual_primary_dist_hist.GetNbinsX();i++){
		double x=manual_primary_dist_hist.GetXaxis()->GetBinCenter(i);
		double y=data->Eval(x);
		pass->SetPoint(pass->GetN(),x,y);
	}
	set_primary_dist(pass);
	delete pass;
}}

void exp_core::set_ruthmask(double lower,double upper,bool reco_v,bool ejec_v){
	TFormula* f=new TFormula("ft","1/(sin(x/2)*sin(x/2)*sin(x/2)*sin(x/2))");
	this->set_primary_dist(f);
	delete f;
	delete this->mask_manual(lower,upper,reco_v,ejec_v);
}

double* exp_core::mask_manual(double lower,double upper,bool reco_v,bool ejec_v){
	double* r = this->generate_mask(lower,upper,ejec_v,reco_v);
	this->apply_mask();
	return r;
}

void exp_core::auto_rutherford(int reps,bool obstructions){
	cout<<endl<<"Accurate up to CM scattering angle "<<happy_ruth_theta(beam_A,beam_Z,targ_A,targ_Z,KE_0_tot_CoM)<<" "<<flush;
	double ranges[14]={0.000001,0.00001,0.0001,0.001,0.01,0.1,0.2,0.3,0.4,0.6,1.0,1.5,2.5,pi};
	vector< double > ejectilecross;
	vector< double > recoilcross;
	for(int i=0;(unsigned)i<detectors.size();i++){
		ejectilecross.push_back(0.0); 
		recoilcross.push_back(0.0);

		this->set_valid_dets(0,i);
		this->set_valid_dets(1,i);
	}

	for(int i=1;i<14;i++){
		vector< int > ejectilecout;
		vector< int > recoilout;
		for(int j=0;(unsigned)j<detectors.size();j++){
			ejectilecout.push_back(0);
			recoilout.push_back(0);
		}
		
		double crosssecforbit=set_rutherford(ranges[i-1],ranges[i]);
// 		this->SetPrimaryDist(0,ranges[i-1],ranges[i]);
		
		for(int j=0;j<reps;j++){
			this->gen_event();
			this->det_check_hits_all(obstructions);
		
			for(int k=0;(unsigned)k<det_hits[0].size();k++){
				recoilout[det_hits[0][k]]++;
			}
		
			for(int k=0;(unsigned)k<det_hits[1].size();k++){
				ejectilecout[det_hits[1][k]]++;
			}
		}

		for(int j=0;(unsigned)j<detectors.size();j++){	
		      recoilcross[j]+=crosssecforbit*(double)recoilout[j]/(double)reps;	
			ejectilecross[j]+=crosssecforbit*(double)ejectilecout[j]/(double)reps;
		}
	}


	double pnA_mb=targ.number_density()*6.24150934E9;	
	cout<<endl<<endl<<setw(14)<<" "<<setw(14)<<"Recoil (mb)"<<setw(14)<<"Ejectile (mb)";
	if(targ.GetThickness()>0)cout<<setw(16)<<"Recoil (Hz/pnA)"<<setw(18)<<"Ejectile (Hz/pnA)";
	for(int i=0;(unsigned)i<detectors.size();i++){
		cout<<endl<<setw(14)<<GetDetName(i)<<setw(14)<<recoilcross[i]<<setw(14)<<ejectilecross[i];
		if(targ.GetThickness()>0)cout<<setw(16)<<recoilcross[i]*pnA_mb<<setw(18)<<ejectilecross[i]*pnA_mb;
		
	}cout<<endl;	
}

void exp_core::auto_rutherford_OTT(int reps,bool obstructions){

	this->set_target_interaction(3);
	cout<<endl<<beam_A<<nuclear_data_ob::get_symb(beam_Z)<<" beam rutherford scattering on "<<targ_A<<nuclear_data_ob::get_symb(targ_Z)<<" target nuclei.";
	this->auto_rutherford(reps,obstructions);cout<<endl;
	
	if(targ.GetThickness()>0&&targ.GetBackThickness()>0){
		TVector3 stopg=targ.GetNormal().Unit()*-stopper_seperation;
		
		offset_primary+=stopg;
		this->set_targ(targ.inverse());
		this->set_target_interaction(3);	
		cout<<endl<<"BACKING SCATTERING";
		cout<<endl<<beam_A<<nuclear_data_ob::get_symb(beam_Z)<<" beam rutherford scattering on "<<targ_A<<nuclear_data_ob::get_symb(targ_Z)<<" backing nuclei.";
		this->auto_rutherford(reps,obstructions);cout<<endl;
		this->set_targ(targ.inverse());
		offset_primary-=stopg;
	}
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void exp_core::gen_uniform(double& theta,double& phi){
	//Generate a random phi theta trajectory evenly distributed on a sphere
	phi=pi*rand.Uniform(0,2);
	theta=acos(rand.Uniform(simm_dist_parB,simm_dist_parA));
}


void exp_core::gen_ruth(double& theta,double& phi){
	//Generate a rutherford phi theta trajectory over the set range
	phi=pi*rand.Uniform(0,2);
	double rr=rand.Uniform(simm_dist_parB,simm_dist_parA);
	theta=2*asin(1/sqrt(rr));
}

void exp_core::gen_manual_primary(double& theta,double& phi){
	phi=pi*rand.Uniform(0,2);
	theta=manual_primary_dist_hist.GetRandom();
}


void exp_core::gen_zeroid(double& theta,double& phi){
	//Generate a random phi theta trajectory evenly distributed on a sphere
	phi=pi*rand.Uniform(0,2);
	theta=2*pi;
	while(theta>pi){theta=abs(rand.Gaus(0,simm_dist_parA));}
}

void exp_core::reset_mask(){
	for(int i=1;i<=manual_primary_dist_mask.GetNbinsX();i++){
		manual_primary_dist_mask.SetBinContent(i,1);
	}
}

void exp_core::apply_mask(){
	//Apply the mask which is all 1 unless it has been manually set
	
	for(int i=1;i<=manual_primary_dist_mask.GetNbinsX();i++){
		int j=i;
		if(reverse_primary)j=(manual_primary_dist_mask.GetNbinsX()+1)-i;
		double y =manual_primary_dist_store.GetBinContent(i)*manual_primary_dist_mask.GetBinContent(j);
		manual_primary_dist_hist.SetBinContent(i,y);
	}
	manual_primary_dist_hist.ComputeIntegral();
}

double* exp_core::generate_mask(double lower,double upper,bool ejec_v,bool reco_v){
	higher_jd(lower,upper);
	if(lower<0)lower=0;
	if(upper>pi)upper=pi;

	double m[2]={ejec_mass,reco_mass+(reco_E_star/jam_phys_amu)};
	bool b[2]={ejec_v,reco_v};
	double u[4]={-1,-1,-1,-1};
	double l[4]={-1,-1,-1,-1};
	
	double P_CoM=momentum_energysplit_CoM(this->get_KE_1_tot_CoM(),m[0],m[1]);
	double betacm=this->get_beta_CoM();
	
	
	//note for some reason I have done ejectile first despite the fact it is second elsewhere in the code
	for(int i=0;i<2;i++){//for the two kinds of particle
		if(b[i]){
			//loop to determine the possible lab angles
			double* ret=lab_boost_CMP_query(betacm,lower,P_CoM,m[i]);
			if(ret[4]>lower){
				//because the distribution is defined for ejectile, there is a pi-x reversal for recoil
				l[i*2]=abs(pi*i-ret[0]);l[i*2+1]=abs(pi*i-ret[2]);
				if(ret[4]<upper){
					u[i*2]=l[i*2+1];u[i*2+1]=l[i*2];
				}else{
					delete ret;ret =lab_boost_CMP_query(betacm,upper,P_CoM,m[i]);
					u[i*2]=abs(pi*i-ret[0]);u[i*2+1]=abs(pi*i-ret[2]);
				}
			}
			delete ret;
		}
	}
	
	double* r=new double[8];
	for(int i=0;i<4;i++){
		if(l[i]>u[i])swap_jd(u[i],l[i]);
		if(!(u[i]>0&&l[i]>0)){
			u[i]=-1;
			l[i]=-1;
		}	
		r[i]=l[i];
		r[i+4]=u[i];
	}
	
	this->reset_mask();
	
	//I considered implementing this with a TF1 but the reversal option for excited beams was a hassle 
	for(int i=1;i<=manual_primary_dist_mask.GetNbinsX();i++){
		bool valid=false;
		double U=manual_primary_dist_mask.GetXaxis()->GetBinUpEdge(i);
		double L=manual_primary_dist_mask.GetXaxis()->GetBinLowEdge(i);
		for(int j=0;j<4;j++){
			double r=abs(u[j]-l[j]);
			if((abs(U-u[j])<r&&abs(U-l[j])<r)||(abs(L-u[j])<r&&abs(L-l[j])<r))valid=true;			
		}
		if(!valid)manual_primary_dist_mask.SetBinContent(i,0);
	}
	
		
	// This is a little more involved because we dont want to just add the two together as that may double count xsec
	// We need to check for overlap in their CM angular ranges
	
	for(int i=0;i<4;i++){
		for(int j=i+1;j<4;j++){// compare each pair to all above it
			if((r[i+4]>r[j])&&(r[i]<r[j+4])){// If there is any overlap
				higher_jd(r[i+4],r[j+4]);// copy the outer limits to one entry
				higher_jd(r[j],r[i]);
				r[i]=-1;// delete the other
				r[i+4]=-1;
				j=4;
			}
		}		
	}
	// Now we have dealt with any overlaps
	
	// It is possible any remaining separate ranges are not ordered
	// This is a problem for one of the following drawing commands	
	
	return r;
}
