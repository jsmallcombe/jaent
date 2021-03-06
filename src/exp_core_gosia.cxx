
//
//
// James's not Geant4 : jaent V1.3.1
// james.smallcombe@outlook.com 01/4/2015
//
//


#include "exp_core.h"

void DecimalOut(ofstream& out,double number){
	if(number==(int)number){
		out<<(int)number<<". ";
	}else{
		out<<number<<" ";
	}
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// Public Members Functions  ///////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

//some default cases
void exp_core::BuildOPINTI(int Nmesh,int NIntegrate,double thetamax){
	
	bool SimpleMesh=false;
	if(Nmesh<0){
		SimpleMesh=true;
		Nmesh=-Nmesh;
		// Introduced this routine because too many points caused issue real gamma detector calcs
		// The phi mesh points should have minimal effect if gamma detectors are phi symmetric anyway
	}
	
	
	thetamax*=pi/180.;
	
	NIntegrate+=NIntegrate%2;
	if(NIntegrate>100)NIntegrate=100;
	if(Nmesh>20){
		cout<<endl<<" Nmesh>20 IS OVERKILL, INCREASE NIntegrate FREELY"<<endl;
	}
	
	double TMin=0;
	double TMax=pi;
	FindThetaMinMax(TMin,TMax);
	if(thetamax>0&&thetamax<TMax)TMax=thetamax;
	double dT=TMax-TMin;

	//To avoid entering a phi range of 0 the theta meshpoints are actually moved in a tiny percent from min/max
	TMin+=dT*0.001;
	TMax-=dT*0.001;


	vector<double> MeshThetaList;
	vector<vector<double>> MeshPhiList;
		
	ShapeToPhiList(MeshThetaList,MeshPhiList,Nmesh,TMin,TMax);
	
	
// 	double Emin=targ.;
	double Emax=targ.beam_e_centre(beam_Z,beam_A,E_beam,TVector3(0.0,0.0,1.0),0.0);
// 	double Emax=dist_target_KE_beam.Eval(0.0);//Needs thick target to have been set
	double Emin=targ.beam_e_centre(beam_Z,beam_A,E_beam,TVector3(0.0,0.0,1.0),1.0);
// 	double Emin=dist_target_KE_beam.Eval(1.0);//Needs thick target to have been set
	int NE=5;
	
	double Emid=(Emin+Emax)*0.5;
	
	if(NE>20)NE=20;
	
	ofstream outI("OPINTI.txt");
	outI<<"OP,INTI"<<endl<<NE<<" "<<-Nmesh<<" "<<Emin<<" "<<Emax<<" "<<TMin*180./TMath::Pi()<<" "<<(TMax)*180./TMath::Pi()<<endl;
	
	vector<double> Epoints;
	for(int e=0;e<NE;e++){
		double E=floor(Emin)+e*(ceil(Emax)-floor(Emin))/(NE-1);
		Epoints.push_back(E);
		DecimalOut(outI,E);
	}
	outI<<endl;	
	
	for(unsigned int i=0;i<MeshThetaList.size();i++){
		DecimalOut(outI,MeshThetaList[i]);
	}
	outI<<endl;
	for(unsigned int i=0;i<MeshPhiList.size();i++){
		int pairs=MeshPhiList[i].size()/2;
		
		if((SimpleMesh&&pairs>2)||pairs>4){
			MeshPhiList[i]=vector<double>{-180,180};
		}
		outI<<pairs<<endl;
		for(unsigned int j=0;j<MeshPhiList[i].size();j++){
			DecimalOut(outI,MeshPhiList[i][j]);
		}
		outI<<endl;
	}
	
	vector<double> IntegThetaList;
	vector<vector<double>> IntegPhiList;
	ShapeToPhiList(IntegThetaList,IntegPhiList,NIntegrate+1,TMin,TMax);
	//NIntegrate+1 because N is segments, so N+1 boundary points

	outI<<NE<<endl;
	for(int e=0;e<NE;e++){DecimalOut(outI,Epoints[e]);}
	outI<<endl;
	for(int e=0;e<NE;e++){outI<<targ.dedx(beam_Z,beam_A,Epoints[e])<<" ";}
	outI<<endl<<NE*2<<" "<<-NIntegrate<<endl;
	
	double dt=(IntegThetaList[1]-IntegThetaList[0])*0.5;
	double meansum=0,phisum=0,meansumruth=0,phisumruth=0;
	
	for(unsigned int i=0;i<IntegPhiList.size();i++){
		double sum=0;
		for(unsigned int j=0;j<IntegPhiList[i].size()-1;j+=2){
			sum+=IntegPhiList[i][j+1]-IntegPhiList[i][j];
		}
		
		double thetamin=IntegThetaList[i]-dt;
		double thetamax=IntegThetaList[i]+dt;
		double scale=1.0;
		
		if(i==0){
			thetamin=IntegThetaList[i];
			scale=0.5;
		}else if(i==IntegPhiList.size()-1){
			thetamax=IntegThetaList[i];
			scale=0.5;
		}
		
		phisum+=sum*scale;
		meansum+=sum*IntegThetaList[i]*scale;
		
		thetamin*=TMath::Pi()/180.0;
		thetamax*=TMath::Pi()/180.0;
		
		double crossthetaterm=(2/pow(sin(0.5*thetamin),2))-(2/pow(sin(0.5*thetamax),2));
		phisumruth+=sum*crossthetaterm;
		meansumruth+=sum*IntegThetaList[i]*crossthetaterm;
	
		DecimalOut(outI,sum);
	}
	
	
	meansum/=phisum;
	meansumruth/=phisumruth;
	
	outI<<endl<<endl<<"Mean Theta: "<<meansum<<endl<<"Rutherford Weighted Mean Theta: "<<meansumruth<<endl;
	outI<<"Rutherford Cross : "<<rutherford_crosssection_lab(beam_mass,targ_mass,targ_Z,beam_Z,Emid,TMin,TMax,0)<<endl;
	
	meansum*=pi/180;
	meansumruth*=pi/180;
	
	double* res=kinetic_lab_calcs_E(Emid,beam_mass,targ_mass,meansum,targ_mass,beam_mass);
	outI<<endl<<"Target Detection"<<endl<<"Mean Theta: "<<-res[5]*180/pi;
	if(res[11]>0)outI<<" OR "<<-res[11]*180/pi;
		
	res=kinetic_lab_calcs_E(Emid,beam_mass,targ_mass,meansumruth,targ_mass,beam_mass);
	outI<<endl<<"Rutherford Weighted Mean Theta: "<<-res[5]*180/pi;
	if(res[11]>0)outI<<" OR "<<-res[11]*180/pi;
	
	outI<<endl<<"Rutherford Cross : "<<rutherford_crosssection_lab(beam_mass,targ_mass,targ_Z,beam_Z,Emid,TMin,TMax,1)<<flush;


// 	kinetic_lab_calcs_readout();
	
	double deltathetaphi=0;
	for(int i=0;(unsigned)i<detectors.size();i++){
		for(int p=0;p<detectors[i].Ntpp();p++){
			deltathetaphi+=detectors[i].tpp(p)->Integral();
		}
	}
	cout<<endl<<" deltathetaphi "<<deltathetaphi<<endl;
	
	outI<<endl;
	outI.close();
}

		
	
void exp_core::ReadDrawMesh(string filename,double Z){
	
	
	ifstream in(filename);
	if(!in.is_open())return;
	
	TGraph* G=new TGraph();
	TGraph* GT=new TGraph();
	
	vector<double> thetav;
	vector<int> Nv;
	vector<double> phi;
	stringstream S;
	string s;
	
	std::getline(in,s);
	if(s.find("INTI")<s.size()){
		std::getline(in,s);
	}
	std::getline(in,s);

	std::getline(in,s);
	S<<s;
	
	double T;
	while(S>>T){
		thetav.push_back(T);
// 		cout<<T<<endl;
	}
	
	int N=0;
	bool inN=true;
	unsigned int count=0;
	while(std::getline(in,s)){
		stringstream SS;
		SS<<s;
		if(inN){
			SS>>N;
			Nv.push_back(N);
		}else{
			double x;
			while(SS>>x){
				phi.push_back(x);
			}
			count++;
		}
		inN=!inN;
		
		if(count==thetav.size())break;
	}
	
	if(thetav.size()!=Nv.size()){
		cout<<endl<<"MISMATCH";
		return;
	}

	unsigned int pn=0;
	for(unsigned int i=0;i<thetav.size();i++){
		double r=tan(thetav[i]*TMath::Pi()/180.)*Z;
		
		for(int n=0;n<Nv[i];n++){
			if(pn<phi.size()-1){
				double a=phi[pn];pn++;
				double b=phi[pn];pn++;
				
				for(double p=a;p<b;p+=0.1){
					double P=p*TMath::Pi()/180.0;
					double x=r*cos(P);
					double y=r*sin(P);
					G->SetPoint(G->GetN(),x,y);
					GT->SetPoint(GT->GetN(),thetav[i],p);
				}
			}
		}
	}

	G->SetMarkerStyle(20);
	GT->SetMarkerStyle(20);
	
	TCanvas* C1 =new TCanvas("A","A",800,800);
	TCanvas* C2 =new TCanvas("B","B",800,800);

	G->SetMarkerStyle(20);
	GT->SetMarkerStyle(20);
	
	C1->cd();
	G->Draw("ap");
	G->GetXaxis()->SetTitle("X [mm]");
	G->GetYaxis()->SetTitle("Y [mm]");
	C2->cd();
	GT->Draw("ap");
	GT->GetXaxis()->SetTitle("#theta Theta Lab [deg]");
	GT->GetYaxis()->SetTitle("#phi Phi [deg]");
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////



void exp_core::ShapeToPhiList(vector<double>& ThetaList,vector<vector<double>>& PhiList,int NPoints,double TMin,double TMax){
	double dT=TMax-TMin;
	double Theta=TMin;
	
	for(int n=0;n<NPoints;n++){
		
		bool trackinside=false;
		
		vector<double> pht;
		for(int i=-1800;i<=1800;i++){
			
			double phi=i*pi/1800.;
			
			bool inside=all_hit_manual(Theta,phi);
			
			if(inside!=trackinside){
				trackinside=!trackinside;
				pht.push_back(i/10.);
			}
		}
		if(trackinside)pht.push_back(180.);
		
		if(pht.size()){
			ThetaList.push_back(Theta*180./TMath::Pi());
			PhiList.push_back(pht);
		}
		Theta+=dT/(NPoints-1);
	}
}


void exp_core::FindThetaMinMax(double& tMin,double& tMax){
	if(use_offsets){
		for(int t=0;t<=1800;t++){
			double T=t*pi/1800.;
			for(int p=0;p<=360;p++){
				double P=p*pi/180.;
				if(all_hit_manual(T,P)){
					tMin=T;
					t=1800;
					break;
				}
			}
		}
		
		for(int t=1800;t>=0;t--){
			double T=t*pi/1800.;
			for(int p=0;p<=360;p++){
				double P=p*pi/180.;
				if(all_hit_manual(T,P)){
					tMax=T;
					t=0;
					break;
				}
			}
		}	
	}else{
		vector<double> minmax;
		for(unsigned int det=0;det<detectors.size();det++){
			minmax.push_back(detectors[det].ThetaMax());
			minmax.push_back(detectors[det].ThetaMin());
		}
		std::sort(minmax.begin(),minmax.end());
		tMin=minmax.front();
		tMax=minmax.back();
	}
}
