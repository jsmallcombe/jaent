
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

void exp_core::draw_exp(int select,int fillo){if(gPad){
	gPad->cd();

	this->draw_fill(select);
	if(select>0){
		twodee->Draw();
		for(int i=0;(unsigned)i<detectors.size();i++){
			if(fillo==2){
				int c=valid_dets[0][i]+valid_dets[1][i]+valid_dets[2][i]+valid_dets[3][i];
				if(c==0)c=12;
				detectors[i].ddraw2->SetFillColor(1+c);
			}else{
				detectors[i].ddraw2->SetFillColor(0);
			}
			if(fillo>0)detectors[i].ddraw2->Draw("sameF");	
			detectors[i].ddraw2->Draw("sameLine");	
		}
	}else{
		threedee->Draw();
		this->draw_arrow();
		for(int i=0;(unsigned)i<detectors.size();i++)detectors[i].ddraw3->Draw("same");
		threedee->Draw("same");//Makes zooming work
	}	

	gPad->Update();
}}

void exp_core::draw_phi(bool fill){if(gPad){
	gPad->cd();
	tp_proj->Draw();
	for(int i=0;(unsigned)i<detectors.size();i++){
		if(fill){detectors[i].tp()->SetFillColor(i+2);
			detectors[i].tp()->Draw("sameF");
			
			for(int p=0;p<detectors[i].Ntpp();p++){
				detectors[i].tpp(p)->SetFillStyle(3004);
				detectors[i].tpp(p)->Draw("sameF");
			}
		}
		else detectors[i].tp()->Draw("same");
	}
	gPad->Update();
}}

void exp_core::draw_xz_labels(){if(gPad){
	gPad->cd();
	this->draw_exp(2);
	for(int i=0;(unsigned)i<detectors.size();i++){
		TText* label= new TText();
		label->SetTextAlign(22);
		label->SetTextSize(0.06);	
		stringstream SS;
		SS<<i;
		label->DrawText(-detectors[i].Xn(), -detectors[i].Zn(), SS.str().c_str());

	}
	gPad->Update();
}}

void exp_core::draw_boost_detectors(bool recoil){if(gPad){
	gPad->cd();	
	double p0=P_1_CoM;
	double m0=reco_mass;
	if(!recoil) m0=ejec_mass;
	
	TH3D drawhist("drawhist","drawhist",200,-worldsize,worldsize,200,-worldsize,worldsize,200,-worldsize,worldsize);

	
	//loop over all detectors (obs)
	for(int i=0;(unsigned)i<detectors.size();i++){
		double t,p;
		TVector3 point;
		detector* d=&detectors[i];
		double magg=sqrt(d->X()*d->X()+d->Y()*d->Y()+d->Z()*d->Z());
		for(int j=0;j<detectors[i].tp()->GetN();j++){
			detectors[i].tp()->GetPoint(j,t,p);
			double* ret=lab_boost_CMP_query(beta_CoM,t,p0,m0);

			//two solutions if in inverse kinematics
			point.SetMagThetaPhi(magg, ret[0], p);
			drawhist.Fill(point.X(),point.Z(),point.Y());
			point.SetMagThetaPhi(magg, ret[2], p);
			drawhist.Fill(point.X(),point.Z(),point.Y());			
			delete ret;
		}		
	}
	drawhist.DrawCopy();
	this->draw_arrow();
	gPad->Update();
}}


void exp_core::draw_target_interaction(int mult,bool obstructions,int reps){if(gPad){
	gPad->cd();	
	
	TH1D interaction("interaction","interaction",302,-1.01,2.02);
	TH1D recoil_decay("recoil_decay","recoil_decay",302,-1.01,2.02);
	recoil_decay.SetLineColor(2);
	interaction.GetXaxis()->SetTitle("Target Fraction");
	
	for(int i=0;i<reps;i++){	
		this->gen_event();
		if((mult>0&&mult<5)||(obstructions&&mult==0)){
			this->det_check_hits_all(obstructions);
			if(this->current_multiplicity()>=mult){
				interaction.Fill(current_target_fraction);
				recoil_decay.Fill(decay_target_fraction);
				i+=100;
			}
		}else{
			interaction.Fill(current_target_fraction);
			recoil_decay.Fill(decay_target_fraction);
		}	
	}

	gPad->SetLogy();interaction.SetMinimum(10);
	interaction.DrawCopy();
	recoil_decay.DrawCopy("same");
	gPad->Update();
}}

void exp_core::draw_decay_Z(int mult,bool obstructions,double zoom){if(gPad&&decay_events){
	gPad->cd();	
	
	if(zoom<=0)zoom=worldsize;
	
	TH1D decay_Z("decay_Z","decay_Z",4000,-zoom,zoom);
	decay_Z.SetLineColor(2);
	decay_Z.GetXaxis()->SetTitle("World Z [mm]");
	
	for(int i=0;i<10000000;i++){	
		this->gen_event();
		if((mult>0&&mult<5)||(obstructions&&mult==0)){
			this->det_check_hits_all(obstructions);
			if(this->current_multiplicity()>=mult){
				decay_Z.Fill(offset_decay.Z());
				i+=100;
			}
		}else{
			decay_Z.Fill(offset_decay.Z());
		}	
	}

	gPad->SetLogy();decay_Z.SetMinimum(10);
	decay_Z.DrawCopy();
	gPad->Update();
}}

TH1D exp_core::theta_cover_auto_norm(bool obstructions,int reps){
	TH1D thetaauto=theta_cover_auto(obstructions,reps);
	for(int i=1;i<=thetaauto.GetNbinsX();i++){
		double N=thetaauto.GetBinContent(i);
		double x1=thetaauto.GetXaxis()->GetBinLowEdge(i);
		double x2=thetaauto.GetXaxis()->GetBinLowEdge(i+1);
		double solid=2*pi*abs(cos(x1)-cos(x2));
		solid/=4*pi;
		thetaauto.SetBinContent(i,N/solid);
	}
	if(gPad){gPad->cd();thetaauto.DrawCopy();gPad->Update();}
	return thetaauto;
}
	
TH1D exp_core::theta_cover_auto(bool obstructions,int reps){
	gPad->cd();
	
	this->set_gun();
	reverse_primary=false;
	this->set_uniform(0,pi);
	this->decay_off();
	
	TH1D thetaauto=theta_hist(reps,obstructions,0,1,0,0);	
	thetaauto.SetName("Theta_coverage");
	thetaauto.SetTitle("Theta_coverage");
	
	return thetaauto;
}

TH1D exp_core::theta_hist(int reps,bool obstructions,bool a,bool b,bool c,bool d){
	bool p[4]={a,b,c,d};
	
	bool hit_check=true;
	if(a+b+c+d<1)hit_check=false;
	
	string s=this->channelnames(ejec_Z,ejec_A);
	if(reverse_primary) s=this->channelnames(reco_Z,reco_A);
	
	s+="_CM_Theta";
	
	TH1D thetah(s.c_str(),s.c_str(),720,0,pi);	
	
	for(int i=0;i<reps;i++){
		this->gen_event();	

		if(hit_check){
			this->det_check_hits_all(obstructions);
			int mult=0;
			for(int j=0;j<4;j++)mult+=p[j]*det_hits[j].size();
			if(mult>0){if(reverse_primary)thetah.Fill(pi-thetaM);else thetah.Fill(thetaM);}
		}else{
			if(reverse_primary){thetah.Fill(pi-thetaM);}else{thetah.Fill(thetaM);}
		}
		
	}	
	
	thetah.Scale(1./reps);

	if(gPad){gPad->cd();thetah.DrawCopy();gPad->Update();}
	return thetah;
}


void exp_core::draw_primary_kinematics(){if(gPad){
	gPad->cd();
	TPad* pad =(TPad*) gPad;//pointer to the pad
	
	double beta = this->get_beta_CoM();
	TVector3 p(0,0,this->get_P_1_CoM());
	
	string es=this->channelnames(ejec_Z,ejec_A);
	string rs=this->channelnames(reco_Z,reco_A);	

	TGraph ang1,ang2,angang,energy1,energy2,ee,beta1,beta2,beta3,beta4,anginvang;
	
	for(int i=0;i<1001;i++){
		double a=pi*i/1000.0;
		TLorentzVector rl=make_lorentzvec(p,reco_mass+reco_E_star/jam_phys_amu);
		TLorentzVector el=make_lorentzvec(-p,ejec_mass);
		el.Boost(0,0,beta);
		rl.Boost(0,0,beta);
		
		ang1.SetPoint(ang1.GetN(),a,rl.Theta());
		ang2.SetPoint(ang2.GetN(),pi-a,el.Theta());
		angang.SetPoint(angang.GetN(),rl.Theta(),el.Theta());
		
		energy1.SetPoint(energy1.GetN(),rl.Theta(),get_KE(&rl));
		energy2.SetPoint(energy2.GetN(),el.Theta(),get_KE(&el));

		beta1.SetPoint(beta1.GetN(),rl.Theta(),get_beta_KE(get_KE(&rl),reco_mass+reco_E_star/jam_phys_amu));
		beta2.SetPoint(beta2.GetN(),el.Theta(),get_beta_KE(get_KE(&el),ejec_mass));
		beta3.SetPoint(beta3.GetN(),rl.Theta(),get_beta_KE(get_KE(&el),ejec_mass));
		beta4.SetPoint(beta4.GetN(),el.Theta(),get_beta_KE(get_KE(&rl),reco_mass+reco_E_star/jam_phys_amu));
		
		anginvang.SetPoint(anginvang.GetN(),pi-a,rl.Theta());
		
		ee.SetPoint(ee.GetN(),get_KE(&rl),get_KE(&el));
		
		p.RotateX(pi/1000.0);
	}

	pad->Divide(3,2);
	pad->cd(1);pad->Update();//need because using drawclone
	ang1.SetMaximum(pi);ang1.GetXaxis()->SetLimits(0,pi);
	ang1.GetXaxis()->SetTitle("CM Angle");	ang1.GetYaxis()->SetTitle("Lab Angle");
	ang1.SetTitle((rs+" red. "+es+" blue.").c_str());
	ang1.SetLineWidth(2);ang2.SetLineWidth(2);ang1.SetLineColor(2);ang2.SetLineColor(4);
	ang1.DrawClone("AC");
	ang2.DrawClone("CSAME");
	pad->cd(2);pad->Update();//need because using drawclone
	energy1.SetMaximum(this->get_E_beam());	energy1.SetMinimum(0);energy1.GetXaxis()->SetLimits(0,pi);
	energy1.GetXaxis()->SetTitle("Lab Angle");	energy1.GetYaxis()->SetTitle("Energy [MeV]");
	energy1.SetTitle((rs+" red. "+es+" blue.").c_str());
	energy1.SetLineWidth(2);energy2.SetLineWidth(2);energy1.SetLineColor(2);energy2.SetLineColor(4);
	energy1.DrawClone("AC");
	energy2.DrawClone("CSAME");
	pad->cd(3);pad->Update();//need because using drawclone
	angang.SetLineWidth(2);
	angang.GetXaxis()->SetTitle(("Lab Angle "+rs).c_str());angang.GetYaxis()->SetTitle(("Lab Angle "+es).c_str());
	angang.DrawClone("AC");
	pad->cd(4);pad->Update();//need because using drawclone
	ee.SetLineWidth(2);
	ee.GetXaxis()->SetTitle(("Energy [MeV] "+rs).c_str());ee.GetYaxis()->SetTitle(("Energy [MeV] "+es).c_str());
	ee.DrawClone("AC");
	pad->cd(5);pad->Update();//need because using drawclone
	beta4.SetMaximum(0.2);	beta4.SetMinimum(0);beta4.GetXaxis()->SetLimits(0,pi);
	beta4.GetXaxis()->SetTitle("Lab Angle"); beta4.GetYaxis()->SetTitle("Beta [c]");
	beta4.SetTitle((rs+" red. "+es+" blue.").c_str());
	beta4.SetLineWidth(1);beta4.SetLineColor(5);beta3.SetLineWidth(1);beta3.SetLineColor(5);
	beta4.DrawClone("AC");
	beta3.DrawClone("CSAME");
	beta1.SetLineWidth(2);beta2.SetLineWidth(2);beta1.SetLineColor(2);beta2.SetLineColor(4);
	beta2.DrawClone("CSAME");
	beta1.DrawClone("CSAME");
	pad->cd(6);pad->Update();
	anginvang.SetLineWidth(2);anginvang.SetLineColor(2);
	anginvang.GetXaxis()->SetTitle(("CM Angle "+es).c_str());anginvang.GetYaxis()->SetTitle("Lab Angle");
	anginvang.SetTitle((rs+" red. "+es+" blue.").c_str());
	anginvang.SetMaximum(pi);
	anginvang.DrawClone("AC");
	ang2.DrawClone("CSAME");
	pad->Update();	
}	
}


// Draws all valid hits for all particles specified
void exp_core::draw_hit_pattern_2D(int proj,int reps,bool obstructions,bool multiR,bool multiE,bool multiA,bool multiB){if(gPad){
	
	// C++11CODE //vector< int > parmin={multiR,multiE,multiA,multiB};
	vector< bool > parmin(4,0);parmin[0]=multiR;parmin[1]=multiE;parmin[2]=multiA;parmin[3]=multiB;

	TH2D hits("Hit_Map","Hit_Map",1000,-worldsize,worldsize,1000,-worldsize,worldsize);
							
	for(int i=1;i<=reps;i++){
		this->gen_event();
		this->det_check_hits_all(obstructions);
		
		for(int j=0;j<4;j++){
			if(parmin[j]){
				int M=det_hits[j].size();
				for(int k=0;k<M;k++){
					TVector3 hit=path_vec_calc(j,det_hits[j][k]);
					switch ( proj ){
						case 1:
							hits.Fill(hit.Y(),hit.X());
							break;
						case 2:
							hits.Fill(hit.X(),hit.Z());
							break;
						default:
							hits.Fill(hit.X(),hit.Y());
							break;
					}
				}
			}
		}
	}
	
	hits.DrawCopy("colz");
	gPad->Update();
}}



//Sums energy and theta of all activated particle hits, summed over all detectors
void exp_core::auto_E_theta(int reps,bool obstructions,bool multiR,bool multiE,bool multiA,bool multiB){if(gPad){
	
	// C++11CODE //vector< int > parmin={multiR,multiE,multiA,multiB};
	vector< bool > parmin(4,0);parmin[0]=multiR;parmin[1]=multiE;parmin[2]=multiA;parmin[3]=multiB;

	TH2D hits("E_Theta","E_Theta",720,0,pi,1000,0,500);
							
	for(int i=1;i<=reps;i++){
		this->gen_event();
		this->det_check_hits_all(obstructions);
		
		for(int j=0;j<4;j++){
			if(parmin[j]){
				if(det_hits[j].size()>0){
					TLorentzVector* four=lor_point[j];
					double KE=four->E()-four->M();
					if(KE>0){
						hits.Fill(four->Theta(),KE);						
					}
				}
			}
		}
	}
	
	hits.DrawCopy("colz");
	gPad->Update();
}}



void exp_core::draw_hits_3D(int projection,int multin,bool obstructions,double display_time,int refresh_rate,bool counterpart_misses,double run_time){if(gPad){//only starts if there is an active pad
	TStopwatch tstop;	
	
	int det_hits_only=0;
	if(counterpart_misses)det_hits_only=-1;//set the -1 flag that draws misses if requested
	   
	gStyle->SetOptStat(0);
	
	int lim=2;//save on some unneeded drawing
	if(decay_events)lim=4;
	
	this->kill_lines();//reset lines
	
	this->draw_exp(projection);	
	
	for(int i=0;i<lim;i++){
		if(projection==0){//if 3D 
			part_pl3[i]->Draw("same");//Draw the particle overlay
		}else{
			part_pl[i]->Draw("same");//Draw the particle overlay
		}		
		if(targ_fusion&&i==0){i++;}//save on some unneeded drawing
	}

	int refresh_count=0;
	double timebzero=tstop.RealTime();tstop.Continue();			
	double T=timebzero;
	
	int counter=0;//Main loop, will exit if no valid draw events for 10000000 events generated
	while(counter<10000000){counter++;
		
		this->gen_event();//gen event and check hits (with or without obstructions)
		this->det_check_hits_all(obstructions);
		
		int mult=this->current_multiplicity();
		
		//stop drawing many blanks if everything stopped in target
		if(multin==0&&mult==0){mult--;for(int i=0;i<lim;i++)if(lor_point[i]->P()>0)mult++;}
		
		//Draw loop
		if(mult>multin-1){//if current event at least required multiplicity (including 0) do a draw
			
			this->kill_lines();//reset lines
			refresh_count++;//count 1 draw;
			counter=0;//good to contiune
			
			for(int i=0;i<lim;i++){//particle loop
				
				// Dont bother draw if P=0
				if(lor_point[i]->P()>0){
					// event if particle has zero hits, possible draw (to box edge)
					if((signed)det_hits[i].size()>det_hits_only){		
						this->draw_path_3D(i,projection);//add either a detector hit or a box edge line to the particle histogram
					}				
				}
				
				if(targ_fusion&&i==0){i++;}//save on some unneeded drawing
			}
			
			//update the draws
			if(refresh_count==1){//On screen update
				
				if(display_time<0&&T>timebzero){//if hold turned on (but not first)
					tstop.Stop();
					cout << "Press ENTER to continue...";
					cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
				}
				T=tstop.RealTime();tstop.Continue();// Start counting display time
			}
			
			//(actually) update the draws
			gPad->Modified();
			gPad->Update();
		}
		
		if(refresh_count>=refresh_rate&&refresh_count>=1){//if not enough draws just loop
			while(tstop.RealTime()-T<display_time){tstop.Continue();}//if display time too small, just wait
			
			//once lap time reached
			if(T-timebzero>run_time){//if its the end, then end
				counter=100000000;
			}else{
				this->kill_lines();//reset lines
				refresh_count=0;					
			}
			tstop.Continue();
		}
	}
	gPad->Update();
}}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////


void exp_core::draw_fill(int select)
{	
	if(select>0) {
		switch ( select ){
			case 1:
				twodee->GetXaxis()->SetTitle("Y-axis [mm]");
				twodee->GetYaxis()->SetTitle("Z-axis [mm]");
				break;
			case 2:
				twodee->GetXaxis()->SetTitle("X-axis [mm]");
				twodee->GetYaxis()->SetTitle("Z-axis [mm]");
				break;
			case 3:
				twodee->GetXaxis()->SetTitle("X-axis [mm]");
				twodee->GetYaxis()->SetTitle("Y-axis [mm]");
				break;
			default:twodee->Reset();
		}		


		for(int i=0;(unsigned)i<detectors.size();i++){
			detectors[i].add_draw(select);
		}
	}
}

void exp_core::draw_arrow(){if(gPad){
	TPolyLine3D pl3(2);
	pl3.SetPoint(0,0,0,0);
	pl3.SetPoint(1,0,0,0);
	Float_t p[6];	
	
	p[0]=offsets[0]->X();
	p[1]=worldsize;
	p[2]=offsets[0]->Y();
	p[3]=p[0];
	p[4]=-p[1];
	p[5]=p[2];
	
	pl3.DrawPolyLine(2,p);
	
	TVector3 be= TVector3(0.0,0.0,1.0);
	if(targ.GetThickness()>0){
		if(targ.GetFormal().Angle(be)>0.001){
			be=targ.GetFormal();
		}
	}
	TVector3 strutA=be.Orthogonal();
	strutA.SetMag(worldsize*0.05);	
	strutA.Rotate(pi/4, be); // rotation around "be"
	TVector3 strutB=strutA;
	strutB.Rotate(pi/2, be);
	
	p[0]=offsets[0]->X()+strutA.X();
	p[1]=offsets[0]->Z()+strutA.Z();
	p[2]=offsets[0]->Y()+strutA.Y();
	p[3]=offsets[0]->X()-strutA.X();
	p[4]=offsets[0]->Z()-strutA.Z();
	p[5]=offsets[0]->Y()-strutA.Y();
	pl3.DrawPolyLine(2,p);	
	
	p[0]=offsets[0]->X()+strutB.X();
	p[1]=offsets[0]->Z()+strutB.Z();
	p[2]=offsets[0]->Y()+strutB.Y();
	p[3]=offsets[0]->X()-strutB.X();
	p[4]=offsets[0]->Z()-strutB.Z();
	p[5]=offsets[0]->Y()-strutB.Y();
	pl3.DrawPolyLine(2,p);
	
	p[0]=offsets[0]->X();
	p[1]=worldsize;
	p[2]=offsets[0]->Y();
	
	for(int i=0;i<20;i++){
		p[3]=offsets[0]->X()+cos(i*pi/10)*worldsize*0.03;
		p[4]=worldsize*0.92;
		p[5]=offsets[0]->Y()+sin(i*pi/10)*worldsize*0.03;
		pl3.DrawPolyLine(2,p);
	}
	
}}


///give a vector for a detector hit or a not hit
TVector3 exp_core::path_vec_calc(int part,int det){
	TVector3 vect(0,0,0);
	TLorentzVector* labframe=lor_point[part];
	if(labframe->P()>0){
		if(det>-1&&(unsigned)det<detectors.size()){
			vect=detectors[det].hit_pos_vec(lor_point[part],offsets[part]);
		}else{
			double p=labframe->Phi(),t=labframe->Theta();
			double X=0,Y=0,Z=0;
			double Zlim=worldsize-offsets[part]->Z();
			if(t>pi*0.5)Zlim=worldsize+offsets[part]->Z();
			double Xlim=worldsize-offsets[part]->X();
			if(abs(p)>pi*0.5)Xlim=worldsize+offsets[part]->X();
			double Ylim=worldsize-offsets[part]->Y();
			if(p<0)Ylim=worldsize+offsets[part]->Y();
			
			double XY=abs(Zlim*tan(t));
			Y=XY*sin(p);
			X=XY*cos(p);
			if(abs(X)>Xlim){
				Y*=abs(Xlim)/abs(X);
				X=Xlim*X/abs(X);;
			}
			if(abs(Y)>Ylim){
				X*=abs(Ylim)/abs(Y);
				Y=Ylim*Y/abs(Y);
			}

			Z=sqrt((X*X)+(Y*Y));
			if(Z==0)Z=Zlim;
			else Z=Z/tan(t);
			vect=TVector3(X,Y,Z);
// 			vect=labframe->Vect();
// 			if(vect.Mag()>0)vect.SetMag(sqrt((X*X)+(Z*Z)+(Y*Y)));
// 			cout<<XY<<" "<<Z/vect.Z()<<" "<<Y/vect.Y()<<" "<<X/vect.X()<<endl;
		}
	}
	return vect;
}


void exp_core::draw_path_3D(int part,int proj,int det){	
	
	this->max_distance_hit(part);
	TVector3 v;
	if(stop_distance[part]>0){
		v=lor_point[part]->Vect();
		if(v.Mag()>0){
			v.SetMag(stop_distance[part]);
		}
	}else{
		v=this->path_vec_calc(part,det);
	}


	switch ( proj ){
		case 1:
			part_pl[part]->SetPoint(0,offsets[part]->Y(),offsets[part]->Z());
			part_pl[part]->SetPoint(1,offsets[part]->Y()+v.Y(),offsets[part]->Z()+v.Z());
			break;
		case 2:
			part_pl[part]->SetPoint(0,offsets[part]->X(),offsets[part]->Z());
			part_pl[part]->SetPoint(1,offsets[part]->X()+v.X(),offsets[part]->Z()+v.Z());
			break;
		case 3:
			part_pl[part]->SetPoint(0,offsets[part]->X(),offsets[part]->Y());
			part_pl[part]->SetPoint(1,offsets[part]->X()+v.X(),offsets[part]->Y()+v.Y());
			break;
		default:
			//YZ flipped for 3D drawing
			part_pl3[part]->SetPoint(0,offsets[part]->X(),offsets[part]->Z(),offsets[part]->Y());
			part_pl3[part]->SetPoint(1,offsets[part]->X()+v.X(),offsets[part]->Z()+v.Z(),offsets[part]->Y()+v.Y());
	}	
	
}	

void exp_core::kill_lines(){	
	for(int i=0;(unsigned)i<part_pl3.size();i++){
		part_pl3[i]->SetPoint(0,0,0,0);
		part_pl3[i]->SetPoint(1,0,0,0);
		part_pl[i]->SetPoint(0,0,0);
		part_pl[i]->SetPoint(1,0,0);
	}
}
