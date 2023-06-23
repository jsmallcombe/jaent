
//
//
// James's not Geant4 : jaent V1.3.1
// james.smallcombe@outlook.com 09/09/2018
//
//


#include "detector_class.h"

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// Public Members Functions  ///////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

detector::detector() : energy(),energyG(), mass(), hit_pat(),energtTheta(),name(), detector_world_place(), detector_world_rotate(),detector_norm()
{
	detector_tp = new TGraph();
	detector_xy = new TGraph();
	ddraw3 = new TPolyLine3D();
	ddraw2 = new TGraph();	
}

detector::detector(vector<double> x,vector<double> y,TVector3 place,TRotation rot) : name(),detector_world_place(place), detector_world_rotate(rot)
{
	
// 	cout<<endl<<detector_norm.X()<<" "<<detector_norm.Y()<<" "<<detector_norm.Z()<<" "<<flush;
	detector_norm=rot * TVector3(0.0,0.0,-1);
// 	cout<<endl<<detector_norm.X()<<" "<<detector_norm.Y()<<" "<<detector_norm.Z()<<" "<<flush;
	
	if(place.Mag()>0)detector_norm*=abs(place.Mag()*TMath::Cos(detector_norm.Angle(place)));
	else detector_norm*=0.000001;
// 	cout<<endl<<detector_norm.X()<<" "<<detector_norm.Y()<<" "<<detector_norm.Z()<<" "<<flush;
	
	// Incase someone was lazy with the rotation
	// and effectively put the detector in facing away
	if(detector_norm.Angle(place)<0.5*pi)detector_norm=-detector_norm;

	energy=TH1D("energy","energy",100000,0,1000);
	energyG=TH1D("energyG","energyG",8000,0,4000);
	energy.GetXaxis()->SetTitle("Energy [MeV]");
	energyG.GetXaxis()->SetTitle("Energy [keV]");
	
	energtTheta=TH2D("energtTheta","energtTheta",1000,0,500000,180,0,pi/2);
	
	mass=TH1D("mass","mass",1000,0,200);
	
	plane=(detector_norm.Dot(place));

	// The rest of the code calculates the outlines xy and tp
	
	//using root libraries phi =0 at y=0 x=+ we dont want any detectors to cross phi+- so use a trick
	//PATCHED OUT//NO DETECTOR should cross y=- x=0//PATCHED OUT
	//PATCHED OUT//for annular START POINT should be y=- x=0//PATCHED OUT
	int numberofpoints=500;
	
	detector_tp = new TGraph();
	detector_xy = new TGraph();
	ddraw3 = new TPolyLine3D();
	ddraw2 = new TGraph();	
	
	while(y.size()<x.size())y.push_back(0.0);

	//Draw the 2D detector shape from list of points

	//First measures perimeter
	double per=0;
	for(int i=0;(unsigned)i<x.size();i++){
		double X0=x[i];
		double Y0=y[i];
		
		double XL=0;
		double YL=0;
		
		if((unsigned)i==x.size()-1){
			XL=(x[0]-X0);
			YL=(y[0]-Y0);
		}else{
			XL=(x[i+1]-X0);
			YL=(y[i+1]-Y0);
		}
		
		per+=sqrt(XL*XL+YL*YL);
	}

	//divided perimeter into 500 points ish
	double maxgap=per/numberofpoints;

	vector< double > th,ph;
	//go through points creating the detector maps
	for(int i=0;(unsigned)i<x.size();i++){

		double X0=x[i];
		double Y0=y[i];
		
		double XL=0;
		double YL=0;
		
		if((unsigned)i==x.size()-1){
			XL=(x[0]-X0);
			YL=(y[0]-Y0);
		}else{
			XL=(x[i+1]-X0);
			YL=(y[i+1]-Y0);
		}
		
		double LL=sqrt(XL*XL+YL*YL);
		int steps=ceil(LL/maxgap);
		XL/=steps;
		YL/=steps;
		
		for(int j=0;j<steps;j++){	
			//cords in detector frame
			double XX=X0+XL*j;
			double YY=Y0+YL*j;
			
			TVector3 world=TVector3(XX,YY,0);
			detector_xy->SetPoint(detector_xy->GetN(),XX,YY);
			
			world=(rot*world)+place;
			
			th.push_back(world.Theta());//temp store the tp points
			ph.push_back(happy_phi(world));//0-2pi converted
		}	
	}
	
	//count the number of times outile crosses phi=0
	int count=0;
	vector< int > crossings;
	for(int i=0;(unsigned)i<ph.size();i++){// count the number of times the phi=0 axis is crossed
		int j=i+1;
		if((unsigned)j==th.size())j=0;
		
		if(abs(ph[i]-ph[j])>pi){
			if((ph[j]>1.5*pi)||(ph[i]>1.5*pi)||(ph[j]<0.5*pi)||(ph[i]>0.5*pi)){
				count++;
				crossings.push_back(i);
				crossings.push_back(j);
			}
		}
	}

	//its a simple detector so can just save it.
	if(count==0)*detector_tp=TGraph((int)th.size(),&th[0],&ph[0]);
	
	//its a theta=0 detectector
	if(count==1){
		double tA=th[crossings[0]],tB=th[crossings[1]];
		double pA=ph[crossings[0]];
		double p0=2*pi+0.1;
		if(pA<pi)p0=-0.1;
		
		//solution to a beam degrading upstream detector
		double endcap=-0.1;
		if(detector_norm.Z()>0	){endcap=pi+0.1;}
		
		for(int i=0;(unsigned)i<ph.size();i++){
			detector_tp->SetPoint(detector_tp->GetN(),th[i],ph[i]);
			
			//add extrapoints to force correct area in tp space
			if(i==crossings[0]){
				detector_tp->SetPoint(detector_tp->GetN(),tA,p0);
				detector_tp->SetPoint(detector_tp->GetN(),endcap,p0);
				detector_tp->SetPoint(detector_tp->GetN(),endcap,pi*2-p0);
				detector_tp->SetPoint(detector_tp->GetN(),tB,pi*2-p0);
						
			}
		}
	}

	//it crosses phi=0 axis but not theta=0	
	if(count==2){
		double tA=th[crossings[0]],tB=th[crossings[1]];
		double pA=ph[crossings[0]];
		double tC=th[crossings[2]],tD=th[crossings[3]];
// 		double pC=ph[crossings[2]];
		double p0=1;
		if(pA<pi)p0=-1;
		double d=0.05;//turns out this bit to add "volume" is un-needed, but I'm leaving
		if(tA>tC)d=-0.05;
		
		for(int i=0;(unsigned)i<ph.size();i++){
			detector_tp->SetPoint(detector_tp->GetN(),th[i],ph[i]);			
			
			//add a wrap around strip through non-physical angle space to link two halves
			if(i==crossings[0]){
				detector_tp->SetPoint(detector_tp->GetN(),tA,pi*(1+p0*(1.1-d)));
				detector_tp->SetPoint(detector_tp->GetN(),-0.1+d,pi*(1+p0*(1.1-d)));
				detector_tp->SetPoint(detector_tp->GetN(),-0.1+d,pi*(1-p0*(1.1-d)));
				detector_tp->SetPoint(detector_tp->GetN(),tB,pi*(1-p0*(1.1-d)));		
			}
			
			if(i==crossings[2]){
				detector_tp->SetPoint(detector_tp->GetN(),tC,pi*(1-p0*(1.1+d)));
				detector_tp->SetPoint(detector_tp->GetN(),-0.1-d,pi*(1-p0*(1.1+d)));
				detector_tp->SetPoint(detector_tp->GetN(),-0.1-d,pi*(1+p0*(1.1+d)));
				detector_tp->SetPoint(detector_tp->GetN(),tD,pi*(1+p0*(1.1+d)));		
			}					
		}
	}

		
	//some kind of crazy irregular shape
	if(count>2){
		*detector_tp=TGraph((int)th.size(),&th[0],&ph[0]);
		cout<<endl<<endl<<" DETECTOR CREATION ERROR "<<count<<" PHI AXIS CROSSINGS."<<endl<<endl;
	}
	
	hit_pat=TH2D("hit_pat","hit_pat",200,detector_xy->GetXaxis()->GetXmin(),detector_xy->GetXaxis()->GetXmax(),200,detector_xy->GetYaxis()->GetXmin(),detector_xy->GetYaxis()->GetXmax());
	
	this->set_draw3();
}

detector::detector( const detector &obj){
	
	detector_tp=(TGraph*)obj.detector_tp->Clone();
	detector_xy=(TGraph*)obj.detector_xy->Clone();
	ddraw3=(TPolyLine3D*)obj.ddraw3->Clone();
	ddraw2=(TGraph*)obj.ddraw2->Clone();
// 	obj.detector_xy->Copy(*detector_xy);
	
	detector_world_place=obj.detector_world_place;
	detector_world_rotate=obj.detector_world_rotate;
	detector_norm=obj.detector_norm;
	
	energy=obj.energy;
	energyG=obj.energyG;
	mass=obj.mass;
	hit_pat=obj.hit_pat;
	energtTheta=obj.energtTheta;
	plane=obj.plane;
	name=obj.name;
}

detector& detector::operator=( const detector &obj){
	if(this!=&obj){//to prevent self-assignment errors
		delete detector_tp;delete detector_xy;
		
		detector_tp=(TGraph*)obj.detector_tp->Clone();
		detector_xy=(TGraph*)obj.detector_xy->Clone();
// 		obj.detector_xy->Copy(*detector_xy);
		
		detector_world_place=obj.detector_world_place;
		detector_world_rotate=obj.detector_world_rotate;
		detector_norm=obj.detector_norm;
		
		energy=obj.energy;mass=obj.mass;hit_pat=obj.hit_pat;energtTheta=obj.energtTheta;plane=obj.plane;name=obj.name;
	}
	return (*this); // for cascading assignment
}

detector::~detector()
{
	delete detector_tp;
	delete detector_xy;
	delete ddraw3;
	delete ddraw2;
}

TVector3 detector::hit_pos_vec(TVector3 motion,TVector3* offset){//not a reference so reused variabl :-S
	if(offset!=NULL){
		if(offset->Mag()>0){
			if(detector_world_place.Mag()==0){
				
				
			}
					
			
			double nume=plane-detector_norm.X()*offset->X()-detector_norm.Y()*offset->Y()-detector_norm.Z()*offset->Z();	
			double deno=motion.X()*detector_norm.X()+motion.Y()*detector_norm.Y()+motion.Z()*detector_norm.Z();
			if(nume!=0&&deno!=0){
				double t=nume/deno;
				if(t>0)return motion*t;
			}
			return TVector3(0.0,0.0,0);	
		}
	}
	
	double mag=detector_norm.Mag()/TMath::Cos(detector_norm.Angle(-motion));
	motion.SetMag(mag);
	return motion;		
	
	return this->hit_pos_vec(motion);
}

TVector3 detector::hit_pos_vec(TLorentzVector* labframe,TVector3* offset){
	return this->hit_pos_vec(labframe->Vect(),offset);
}

TVector3 detector::hit_pos_vec(double theta, double phi){
	TVector3 tmp_vec;
	tmp_vec.SetMagThetaPhi(1,theta,phi);
	return this->hit_pos_vec(tmp_vec);
}

bool detector::hit_offset_check(TVector3 motion,TVector3* offset){
	TVector3 relative=this->hit_pos_vec(motion,offset);

	if(relative.Mag()>0){
		TVector3 ondet=relative-detector_world_place;
		if(offset!=NULL)ondet+=(*offset);
		ondet*=detector_world_rotate.Inverse();
			//cout<<relative.Mag()<<endl;
		if(detector_xy->IsInside(ondet.X(),ondet.Y()))return true;
	}
	return false;	
}

bool detector::hit_offset_check(TLorentzVector* labframe,TVector3* offset){
	return this->hit_offset_check(labframe->Vect(),offset);
}

void detector::set_draw3(){
	double x,y;
	ddraw3->SetPolyLine(detector_xy->GetN());
	for(int i=0;i<detector_xy->GetN();i++){
		detector_xy->GetPoint(i,x,y);	
		TVector3 world=TVector3(x,y,0);
		world=(detector_world_rotate*world)+detector_world_place;
		ddraw3->SetPoint(i,world.X(),world.Z(),world.Y());		
	}
}

void detector::add_draw(TH3* addhist){
	double x,y;
	for(int i=0;i<detector_xy->GetN();i++){
		detector_xy->GetPoint(i,x,y);	
		TVector3 world=TVector3(x,y,0);
		world=(detector_world_rotate*world)+detector_world_place;
		addhist->Fill(world.X(),world.Z(),world.Y());		
	}
}


void detector::add_draw(TH2* addhist,int projection){
	double x,y;
	
	for(int i=0;i<detector_xy->GetN();i++){
		detector_xy->GetPoint(i,x,y);	
		TVector3 world=TVector3(x,y,0);
		world=(detector_world_rotate*world)+detector_world_place;
		
		switch ( projection ){
			case 1:addhist->Fill(world.Y(),world.Z());break;
			case 2:addhist->Fill(world.X(),world.Z());break;
			default:addhist->Fill(world.X(),world.Y());
		}				
	}
}

void detector::add_draw(int projection){
	double x,y;
	ddraw2->Clear();
	for(int i=0;i<detector_xy->GetN();i++){
		detector_xy->GetPoint(i,x,y);	
		TVector3 world=TVector3(x,y,0);
		world=(detector_world_rotate*world)+detector_world_place;
		
		switch ( projection ){
			case 1:ddraw2->SetPoint(ddraw2->GetN(),world.Y(),world.Z());break;
			case 2:ddraw2->SetPoint(ddraw2->GetN(),world.X(),world.Z());break;
			default:ddraw2->SetPoint(ddraw2->GetN(),world.X(),world.Y());
		}				
	}
}


void detector::Fill(TLorentzVector* four,TVector3* offset){
	double massu=four->M();
	double KE=four->E()-massu;
	massu/=jam_phys_amu;
	
	TVector3 hit=this->hit_pos_vec(four->Vect(),offset)-detector_world_place;
	if(offset!=NULL)hit+=(*offset);
	hit*=detector_world_rotate.Inverse();
	
	energy.Fill(KE);
	energyG.Fill(KE*1000);
	mass.Fill(massu);
	hit_pat.Fill(hit.X(),hit.Y());
	energtTheta.Fill(KE*1000,four->Vect().Theta());
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////OTHER //////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

TH1D* add_resolution(TH1D* in,double percent,double at){
	string s=in->GetName();
	TH1D* blur=(TH1D*)in->Clone((s+"_blur").c_str());
	blur->Reset();
	TAxis* x=in->GetXaxis();
	TRandom2 rand;rand.SetSeed();
	
	for(int i=1;i<=x->GetNbins();i++){
		double e=in->GetBinCenter(i);
		double de=percent*sqrt(e*at)/100;
		de/=2.355;//convert to sigma from FWHM
		
		for(int j=0;j<in->GetBinContent(i);j++)	blur->Fill(rand.Gaus(e,de));
	}
	
	return blur;
}


