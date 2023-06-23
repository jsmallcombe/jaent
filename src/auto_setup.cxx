
//
//
// James's not Geant4 : jaent V1.3.1
// james.smallcombe@outlook.com 01/4/2015
//
//


#include "auto_setup.h"


void mwpc_auto_setup(exp_core &exp,int select,double silicon_Z,bool silicon_detail){
	exp.reset_detectors();
	
	//define some 2D things
	
	// C++11CODE //vector<double> xa= {-100, -100, 100, 100};
	// C++11CODE //vector<double> ya= {-100, 100, 100, -100};
	// C++11CODE //vector<double> xb= {-200, -200, 200, 200};
	// C++11CODE //vector<double> yb= {-100, 100, 100, -100};
	vector<double> xa;xa.push_back(-100);xa.push_back(-100);xa.push_back(100);xa.push_back(100);
	vector<double> ya;ya.push_back(-100);ya.push_back(100);ya.push_back(100);ya.push_back(-100);
	vector<double> xb;xb.push_back(-200);xb.push_back(-200);xb.push_back(200);xb.push_back(200);
	vector<double> yb;yb.push_back(-100);yb.push_back(100);yb.push_back(100);yb.push_back(-100);
	
	TVector3 position;
	TRotation rotit;	
	
	if(silicon_Z){
		
		add_S1(exp,silicon_Z,silicon_detail);
		
		if(silicon_detail){
			double d_p=pi*(22.5/180)/40;			
			vector<double> xdE,ydE;	
			for(int i=0;i<41;i++){
				double x=sin(pi*(-11.25/180)+((double)i*d_p))*20.9;
				double y=(cos(pi*(-11.25/180)+((double)i*d_p))*20.9)+10.19;
				xdE.push_back(x);ydE.push_back(y);
			}
			for(int i=0;i<41;i++){
				double x=sin(pi*(11.25/180)-((double)i*d_p))*37.3;
				double y=(cos(pi*(11.25/180)-((double)i*d_p))*37.3)+10.19;
				xdE.push_back(x);ydE.push_back(y);
			}	

			//dE
			position=TVector3(0.0,0.0,silicon_Z+20);
			rotit=TRotation();
			//rotit.RotateZ(pi*(22.5/180)/2);
			rotit.RotateX(-0.7365);
			rotit.RotateZ(-pi/12);	

			for(int i=1;i<13;i++){
				if(i!=5&&i!=11)exp.add_detector(xdE,ydE,position,rotit,false,true);
				rotit.RotateZ(-pi/6);
			}
			
			for(int i=2;i<12;i++){
				exp.set_pair(1,i);
			}
			
		}
		
	// 	if(silicon_detail&&silicon_on){
	// 		for(int i=0;i<13;i++){
	// 			exp.set_pair(0,i);
	// 			exp.set_pair(i,0);		
	// 		}
	// 		for(int i=0;i<5;i++){
	// 			exp.set_pair(i+3,1);
	// 			exp.set_pair(1,i+3);
	// 			exp.set_pair(i+8,2);
	// 			exp.set_pair(2,i+8);		
	// 		}
	// 	}
	}
	

	rotit=TRotation();
	switch ( select ){
		case 1:
			position=TVector3(249,0,0);
			rotit.RotateY(TMath::Pi()/2);
			exp.add_detector(xa,ya,position,rotit,false,false,true,true);
			position=TVector3(-249,0,0);
			rotit.RotateY(-TMath::Pi());
			exp.add_detector(xa,ya,position,rotit,false,false,true,true);
			break;
		case 2:
			position=TVector3(215,0,27);
			rotit.RotateY(TMath::Pi()/2);
			exp.add_detector(xa,ya,position,rotit,false,false,true,true);
			position=TVector3(-215,0,27);
			rotit.RotateY(-TMath::Pi());
			exp.add_detector(xa,ya,position,rotit,false,false,true,true);
			break;
		case 3:
			position=TVector3(304.9,0,304.9);
			rotit.RotateY(TMath::Pi()/4);
			exp.add_detector(xb,yb,position,rotit,false,false,true,true);
			position=TVector3(-304.9,0,304.9);
			rotit.RotateY(-TMath::Pi()/2);
			exp.add_detector(xb,yb,position,rotit,false,false,true,true);
			break;
		default:
			position=TVector3(cos(0.873)*224,0,sin(0.873)*224);
			rotit.RotateY(0.873);
			exp.add_detector(xa,ya,position,rotit,false,false,true,true);
			position=TVector3(-cos(0.873)*224,0,sin(0.873)*224);
			rotit.RotateY(-1.746);
			exp.add_detector(xa,ya,position,rotit,false,false,true,true);
			rotit=TRotation();
			position=TVector3(cos(0.873)*224,0,-sin(0.873)*224);
			rotit.RotateY(2.269);
			exp.add_detector(xa,ya,position,rotit,false,false,true,true);
			position=TVector3(-cos(0.873)*224,0,-sin(0.873)*224);
			rotit.RotateY(-4.538);
			exp.add_detector(xa,ya,position,rotit,false,false,true,true);
	}	
}

void add_rutherford_monitor(exp_core &exp,double XX,double ZZ,double RR){
	RR=abs(RR);
	
// 	bool zero=false;
	
	double extraang=0;
	double ang=atan(XX/ZZ);
// 	if(ZZ<0)ang=pi-ang;
	
	/// All this stuff is redundent now the phi=0 code improved in detector class
	if(cos(ang)*RR>abs(XX)){
// 		zero=true;
		double proj=-XX/cos(ang);
		extraang=pi*0.5-acos(proj/RR);
	}
	
	vector<double> x,y;
	for(int i=0;i<360;i++){
	  x.push_back(RR*sin(((double)i*pi/180)+extraang));
	  y.push_back(-RR*cos(((double)i*pi/180)+extraang));
	}	
	
	TRotation rotit;
	rotit.RotateY(ang);
	exp.add_detector(x,y,TVector3(XX,0,ZZ),rotit,true,true);
}

void gen_annulus(vector<double>& x,vector<double>&y,double rr,double RR){
	x.clear();
	y.clear();
	
	RR=abs(RR);
	rr=abs(rr);
	if(rr>RR){double x=rr;rr=RR;RR=x;}
	
	for(int i=0;i<=360;i++){
	  y.push_back(RR*sin(((double)i*pi/180)));
	  x.push_back(-RR*cos(((double)i*pi/180)));
	}	
	for(int i=0;i<=360;i++){
	  y.push_back(-rr*sin(((double)i*pi/180)));
	  x.push_back(-rr*cos(((double)i*pi/180)));
	}	
}
	
void add_annulus(exp_core &exp,double rr,double RR,TVector3 vec,TRotation rot){
	vector<double> x,y;
	gen_annulus(x,y,rr,RR);
	exp.add_detector(x,y,vec,rot,true,true,true);
}
	
void add_annulus(exp_core &exp,double ZZ,double rr,double RR){
	add_annulus(exp,rr,RR,TVector3(0,0,ZZ));
}

void add_micron_S(exp_core &exp,double ZZ,bool detail,int type){
	int N=24;
	double RR=35;	
	double rr=11;
	if(type==1){RR=48;rr=24;N=16;}

	if(detail){
		// C++11CODE //vector<double> xpcb= {0, 40, 60,60,40,-40,-60,-60,-40,0};
		// C++11CODE //vector<double> ypcb= {60,60,40,-40,-60,-60,-40,40,60,60};
		
		vector<double> xpcb; double a[10]={0, 40, 60,60,40,-40,-60,-60,-40,0};
		vector<double> ypcb; double b[10]={60,60,40,-40,-60,-60,-40,40,60,60};
		for(int i=0;i<10;i++){xpcb.push_back(a[i]);ypcb.push_back(b[i]);}
		
		for(int i=200;i>=0;i--){
			double x=RR*sin(pi*((double)i/100));
			double y=RR*cos(pi*((double)i/100));
			xpcb.push_back(x);
			ypcb.push_back(y);
		}		
		double pitch=(RR-rr)/N;	
		for(int i=0;i<N;i++){
			add_annulus(exp,ZZ,rr+(i*pitch)+pitch*0.05,rr+((i+1)*pitch));
			stringstream ss;
			ss<<"S"<<type<<" Ring "<<i;
			exp.SetDetName(ss.str());			
		}
		exp.add_detector(xpcb,ypcb,TVector3(0.0,0.0,ZZ));
		stringstream ss;
		ss<<"S"<<type<<" PCB";
		exp.SetDetName(ss.str());	
        
	}
	else{
		add_annulus(exp,ZZ,rr,RR);
		stringstream ss;
		ss<<"S"<<type<<" Silicon";
		exp.SetDetName(ss.str());
	}
}
void add_S1(exp_core &exp,double ZZ,bool detail){add_micron_S(exp,ZZ,detail,1);}
void add_S2(exp_core &exp,double ZZ,bool detail){add_micron_S(exp,ZZ,detail,2);}
void add_S3(exp_core &exp,double ZZ,bool detail){add_micron_S(exp,ZZ,detail,3);}


void add_PIN(exp_core &exp,TVector3 vec,TRotation rotit,bool centre,bool pcb){
	
	vector<double> xpcb; double a[8]={-5.4,-5.4,5.4,5.4,0.5,0.5,-0.5,-0.5};
	vector<double> ypcb; double b[8]={-7.275,5.125,5.125,-7.275,-7.275,-6.275,-6.275,-7.275};
	for(int i=0;i<8;i++){xpcb.push_back(a[i]);ypcb.push_back(b[i]);}
	
	vector<double> xactive; double aa[4]={-4.865,3.975,3.975,-4.865};
	vector<double> yactive; double bb[4]={4.89,4.89,-4.89,-4.89};
	for(int i=0;i<4;i++){xactive.push_back(aa[i]);yactive.push_back(bb[i]);}
	
	if(!centre){
		for(int i=0;i<8;i++)ypcb[i]+=7.275;
		for(int i=0;i<4;i++)yactive[i]+=7.275;
	}
	
	TVector3 smidge(0,0,-0.5);
	smidge*=rotit;//make detector normal
	if(smidge.Angle(TVector3(0,0,-1))>pi/2)smidge=-smidge;

	if(pcb)exp.add_detector(xpcb,ypcb,vec,rotit);
	exp.add_detector(xactive,yactive,vec+smidge,rotit,true,true,true,true);
	
}

void add_square(exp_core &exp,TVector3 vec,TRotation rotit,double length){
	
	vector<double> x; 
	x.push_back(length/2);x.push_back(length/2);x.push_back(-length/2);x.push_back(-length/2);
	vector<double> y; 
	y.push_back(length/2);y.push_back(-length/2);y.push_back(-length/2);y.push_back(length/2);	

	exp.add_detector(x,y,vec,rotit,true,true,true,true);
	
}

void add_target_ladder(exp_core &exp,double XX,double ZZ,double RR){
	XX/=2;ZZ/=2;
	if(RR>XX)XX=RR+2;
	if(RR>ZZ)ZZ=RR+2;

	// C++11CODE //vector<double> x={XX,XX,-XX,-XX,XX,XX};
	// C++11CODE //vector<double> y={0,ZZ,ZZ,-ZZ,-ZZ,0};
	vector<double> x;x.push_back(XX);x.push_back(XX);x.push_back(-XX);x.push_back(-XX);x.push_back(XX);x.push_back(XX);
	vector<double> y;y.push_back(0);y.push_back(ZZ);y.push_back(ZZ);y.push_back(-ZZ);y.push_back(-ZZ);y.push_back(0);



	for(int i=0;i<360;i++){
		y.push_back(-RR*sin((double)i*pi/180));
		x.push_back(RR*cos((double)i*pi/180));
	}	

	TVector3 place=exp.targ.GetNormal();
	place.SetMag(0.01);
	place+=exp.offset_primary;
	
	TRotation rotit;
	rotit.RotateY(exp.targ.GetNormal().Theta());
	rotit.RotateZ(exp.targ.GetNormal().Phi());
	exp.add_detector(x,y,place,rotit);
}


void add_chamber_cylinder(exp_core &exp,double ZZ,double RR,bool upright){
	ZZ=abs(ZZ);
	RR=abs(RR);
	
	
	vector<double> x,y;
	for(int i=0;i<360;i++){
	  y.push_back(RR*1.16*sin(((double)i*pi/180)));
	  x.push_back(-RR*1.16*cos(((double)i*pi/180)));
	}

	vector<double> X,Y;
	  Y.push_back((ZZ/2)+10);Y.push_back((ZZ/2)+10);Y.push_back(-(ZZ/2)-10);Y.push_back(-(ZZ/2)-10);
	  X.push_back((RR*0.5773)+10);X.push_back(-(RR*0.5773)-10);X.push_back(-(RR*0.5773)-10);X.push_back((RR*0.5773)+10);
	  
	TRotation rotit;
	
	if(upright){
		rotit.RotateX(pi/2);
		exp.add_detector(x,y,TVector3(0,ZZ/2,0),rotit,false,false,false,false);
		exp.add_detector(x,y,TVector3(0,-ZZ/2,0),rotit,false,false,false,false);
		rotit.RotateX(-pi/2);
		for(int i=0;i<6;i++){
			exp.add_detector(X,Y,TVector3(RR*sin(i*(pi/3)),0,RR*cos(i*(pi/3))),rotit,false,false,false,false);
			rotit.RotateY(pi/3);		
		}
	}else{
		exp.add_detector(x,y,TVector3(0,0,ZZ/2),rotit,false,false,false,false);
		exp.add_detector(x,y,TVector3(0,0,-ZZ/2),rotit,false,false,false,false);
		rotit.RotateX(pi/2);
		for(int i=0;i<6;i++){
			exp.add_detector(X,Y,TVector3(RR*sin(i*(pi/3)),RR*cos(i*(pi/3)),0),rotit,false,false,false,false);
			rotit.RotateZ(-pi/3);		
		}
	}

	
}


void add_chamber_cubeoid(exp_core &exp,double ZZ,double YY,double XX){
	ZZ=abs(ZZ);
	YY=abs(YY);
	XX=abs(XX);
	
	vector<double> X,Y;
	Y.push_back((YY/2)+10);Y.push_back((YY/2)+10);Y.push_back(-(YY/2)-10);Y.push_back(-(YY/2)-10);
	X.push_back((XX/2)+10);X.push_back(-(XX/2)-10);X.push_back(-(XX/2)-10);X.push_back((XX/2)+10);
  
	TRotation rotit;
	
	exp.add_detector(X,Y,TVector3(0,0,ZZ/2),rotit,false,false,false,false);
	exp.add_detector(X,Y,TVector3(0,0,-ZZ/2),rotit,false,false,false,false);
	
	Y.clear();
	X.clear();
	
	Y.push_back((ZZ/2)+10);Y.push_back((ZZ/2)+10);Y.push_back(-(ZZ/2)-10);Y.push_back(-(ZZ/2)-10);
	X.push_back((XX/2)+10);X.push_back(-(XX/2)-10);X.push_back(-(XX/2)-10);X.push_back((XX/2)+10);
	
	rotit.RotateX(pi/2);
	
	exp.add_detector(X,Y,TVector3(0,YY/2,0),rotit,false,false,false,false);
	exp.add_detector(X,Y,TVector3(0,-YY/2,0),rotit,false,false,false,false);
	
	Y.clear();
	X.clear();	
	
	Y.push_back((ZZ/2)+10);Y.push_back((ZZ/2)+10);Y.push_back(-(ZZ/2)-10);Y.push_back(-(ZZ/2)-10);
	X.push_back((YY/2)+10);X.push_back(-(YY/2)-10);X.push_back(-(YY/2)-10);X.push_back((YY/2)+10);
	
	rotit.RotateZ(pi/2);
	
	exp.add_detector(X,Y,TVector3(XX/2,0,0),rotit,false,false,false,false);
	exp.add_detector(X,Y,TVector3(-XX/2,0,0),rotit,false,false,false,false);
	
}

void spice_auto_setup(exp_core &exp,int select,int S3,double overr){
	exp.reset_detectors();
	
	double S3_mm=22;
	double PDA_mm=28;
	
	//0: // Raised Target no stopper
	//1: // Raised Target WITH stopper
	//3: // No raised target
	if(select==3){
		S3_mm=18;
		exp.set_target_primary_offset(0,0,0);
		exp.set_stopper_seperation(0);
	}else{
		exp.set_target_primary_offset(0,0,-8);
		if(select==1)exp.set_stopper_seperation(7.5);
	}

	if(select==3){
		add_annulus(exp,0,4,100);
	}else{
		if(select==1)add_annulus(exp,0,8,100);
		else add_annulus(exp,2,9,100);
	}	
	exp.SetDetName("Target Wheel");	
	exp.set_valid_part(0,false,false,false,false);
	
	if(select!=3){
		add_annulus(exp,40,60,TVector3(0,17,11.63));
		exp.set_valid_part(1,false,false,false,false);
		exp.SetDetName("Target Wheel B");	
		
		add_annulus(exp,-8,4,8);
		exp.set_valid_part(2,false,false,false,false);	
		exp.SetDetName("Target Frame");		
	}
	
	if(overr>0){S3_mm=overr;PDA_mm=overr;}
	switch ( S3 ){
		case 1:add_S3(exp,S3_mm,false);break;
		case 2:add_S3(exp,S3_mm,true);break;
		case 3:add_PDA(exp,PDA_mm,false);break;
		default: break;
	}
}

void add_PDA(exp_core& experiment,double Z_nom,bool PCB){
	
	double safety_gap=0.5;
	double A_spacing=10.6+safety_gap;
	double B_spacing=12.4+safety_gap;
	double B_start_offset=12.4+5.3+safety_gap;
	double centre_pixel_shift=6.2;
			
	TRotation rota,rotb;
	rotb.RotateZ(pi);
		
	//Quadrants 
	for(int k=0;k<4;k++){//quad
		double Z=Z_nom;
		double X=0;
		
		for(int i=0;i<4;i++){//axis A
			double Y=B_start_offset;
			
				for(int j=0;j<3;j++){//axis B
					//special rules
					double ys=0;			
					if(i==0){
						ys=centre_pixel_shift;
						if(k%2)ys+=2.0;
					}
					if(i==3 && j==2)continue;
					
					//add detector
					TVector3 vec(-X,Y+ys,Z);
					vec*=rota;
					add_PIN(experiment,vec,rotb,false,PCB);
					
					//step for next row
					Y+=B_spacing;
				}
			//step for next column
			X+=A_spacing;	
		}
		//Rotate for the next quadrant
		rota.RotateZ(pi/2);
		rotb.RotateZ(pi/2);	
	}
}




void AddGapS3Rings(exp_core &exp,TVector3 pos,TRotation rot,int InnerRing,int OuterRing,int OffSectorStart,int OffSectorEnd){

	double phiA=-180.0,phiB=180.0;
	if(OffSectorEnd>=0){
		phiA=(OffSectorEnd+0.5)*11.25;
		phiB=(OffSectorStart-0.5)*11.25+360;
	}
	vector<double> x,y;
	for(int i=phiA*10;i<=phiB*10;i++){
		double ang=i*TMath::Pi()/1800.;
		x.push_back((11.0+InnerRing)*cos(ang));
		y.push_back((11.0+InnerRing)*sin(ang));
	}
	for(int i=phiB*10;i>=phiA*10;i--){
		double ang=i*TMath::Pi()/1800.;
		x.push_back((12.0+OuterRing)*cos(ang));
		y.push_back((12.0+OuterRing)*sin(ang));
	}	
	
	exp.add_detector(x,y,pos,rot,true,true);
}


void AddS3BadPix(exp_core &exp,TVector3 pos,TRotation rot,int InnerRing,int OuterRing,vector<vector<int>> BadPix,int SectorOffsetN){
	
	vector<vector<bool>> OnPix(24,vector<bool>(32,true));
	
	for(int r=0;r<InnerRing;r++){
		for(int s=0;s<32;s++){
			OnPix[r][s]=false;
		}
	}
	
	for(int r=OuterRing+1;r<24;r++){
		for(int s=0;s<32;s++){
			OnPix[r][s]=false;
		}
	}
	
	for(unsigned int i=0;i<BadPix.size();i++){
			int s=BadPix[i][1]+SectorOffsetN;
			if(s>=32)s-=32;
			if(s<0)s+=32;
// 			cout<<endl<<BadPix[i][0]<<" "<<s;
		OnPix[BadPix[i][0]][s]=false;
	}
	
	for(int r=23;r>=0;r--){
		cout<<endl;
		for(int s=0;s<32;s++){
			cout<<OnPix[r][s];
		}
	}	
	
	bool started=false;
	vector<double> x,y;
	int s=0;
	for(s=0;s<32;s++){
		
		bool good=false;
		for(int r=InnerRing;r<=OuterRing;r++){
			if(OnPix[r][s]){
				started=true;
				good=true;
				double ang=(s-0.5)*11.25;
				for(int i=0;i<=40;i++){
					x.push_back((11.0+r)*cos(ang*TMath::Pi()/180.));
					y.push_back((11.0+r)*sin(ang*TMath::Pi()/180.));
					ang+=11.25/40.;
				}
				break;
			}
		}
		
		if(!good&&started)break;
	}
	

	double SGAP=s;
	for(int S=SGAP-1;S>=SGAP-32;S--){
		s=S;
		if(s<0)s+=32;

		bool good=false;
		for(int r=OuterRing;r>=InnerRing;r--){
			if(OnPix[r][s]){
				good=true;
				double ang=(s+0.5)*11.25;					
				for(int i=0;i<=40;i++){
					x.push_back((12.0+r)*cos(ang*TMath::Pi()/180.));
					y.push_back((12.0+r)*sin(ang*TMath::Pi()/180.));
					ang-=11.25/40.;
				}
				break;
			}
		}
		
		if(!good){
			break;
		}
	}
	

	if(s==0){
		exp.add_detector(x,y,pos,rot,true,true);
		return;
	}
	
	double SGAPB=s;
	for(s=SGAPB+1;s<32;s++){
		bool good=false;
		for(int r=InnerRing;r<=OuterRing;r++){
			if(OnPix[r][s]){
				good=true;
				double ang=(s-0.5)*11.25;
				for(int i=0;i<=40;i++){
					x.push_back((11.0+r)*cos(ang*TMath::Pi()/180.));
					y.push_back((11.0+r)*sin(ang*TMath::Pi()/180.));
					ang+=11.25/40.;
				}
				break;
			}
		}
		
		if(!good)break;
	}
	exp.add_detector(x,y,pos,rot,true,true);
	return;
}



