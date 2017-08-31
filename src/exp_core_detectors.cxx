
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

void exp_core::add_detector(vector< vector<double> > points,TVector3 place,TRotation rot,bool A,bool B,bool C,bool D){
	detector det=detector(points[0],points[1],place,rot);
	this->add_detector(det,A,B,C,D);
}

void exp_core::add_detector(vector<double> x,vector<double> y,TVector3 place,TRotation rot,bool A,bool B,bool C,bool D){
	detector det=detector(x,y,place,rot);
	this->add_detector(det,A,B,C,D);
}

void exp_core::add_detector(detector det,bool A,bool B,bool C,bool D)
{
	detectors.push_back(det);

	this->extend_doubles();
	
	valid_dets[0].push_back(A);
	valid_dets[1].push_back(B);
	valid_dets[2].push_back(C);
	valid_dets[3].push_back(D);	
}

void exp_core::set_pair(int a,int b,bool t)
{
	if((unsigned)a<detectors.size()&&(unsigned)b<detectors.size()&&a>-1&&b>-1){
		if(a==b){
			valid_doubles[a][b]=true;
		}else{
			valid_doubles[a][b]=t;
			valid_doubles[b][a]=t;		
		}
	}
}


void exp_core::reset_detectors()
{
	detectors.resize(0);
	threedee->Reset();
	tp_proj->Reset();
	twodee->Reset();	
	valid_dets.assign(4,vector< bool >());
	valid_doubles.resize(0);
}

void exp_core::reset_valid_dets(){valid_dets.assign(4,vector< bool >(detectors.size(),false));}

void exp_core::set_valid_dets(int part,int in){if((unsigned)in<detectors.size()&&part<4)valid_dets[part][in]=true;}

void exp_core::set_valid_part(int in,bool A,bool B,bool C,bool D){
	if((unsigned)in<detectors.size()){
		valid_dets[0][in]=A;
		valid_dets[1][in]=B;
		valid_dets[2][in]=C;
		valid_dets[3][in]=D;
	}
}


void exp_core::reset_doubles()
{
	valid_doubles.resize(0);
	for(int i=0;(unsigned)i<detectors.size();i++)valid_doubles.push_back(vector<bool>(detectors.size(),false));
	for(int i=0;(unsigned)i<detectors.size();i++)valid_doubles[i][i]=true;	
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// PRIVATE Members Functions  //////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////



///////////////////////////
////// Detector hits //////
///////////////////////////

/// check a single detector
bool exp_core::det_hit_quick(int det,int part){
	if((unsigned)det<detectors.size()&&lor_point[part]->P()>0){
		if(use_offsets)if(offsets[part]->Mag()>0){//Only if it fails the first 2 IFs do the second method
			if(detectors[det].hit_offset_check(lor_point[part],offsets[part]))
				return true;
			else 
				return false;
		}
		if(detectors[det].tp()->IsInside(lor_point[part]->Theta(),happy_phi(lor_point[part]->Phi())))return true;
	}
	return false;	
}

/// check all detectors for one particle
void exp_core::det_check_hits(int part,bool obstruct){
	if(part>=0&&part<4){
		det_hits[part].resize(0);
		stop_distance[part]=-1;
		
		if(lor_point[part]->P()>0){
			
			//If no obstruction then just check all the allowed detectors
			if(!obstruct)for(int i=0;(unsigned)i<detectors.size();i++)
				if(valid_dets[part][i])
					if(this->det_hit_quick(i,part))det_hits[part].push_back(i);
			
			//If obstruction check ALL detectors
			if(obstruct){
				vector< int > possible_hits;
				vector< double > possible_distance;
				
				for(int i=0;(unsigned)i<detectors.size();i++){ // create a list of those on path
					if(this->det_hit_quick(i,part)){
						possible_hits.push_back(i);
						possible_distance.push_back(detectors[i].hit_pos_vec(lor_point[part],offsets[part]).Mag());
	// 					if(use_offsets)
	// 					else possible_distance.push_back(detectors[i].hit_pos_vec(lor_point[part]).Mag());
					}
				}
				
				//order the hits based on distance	
				for(int i=1;(unsigned)i<possible_hits.size();i++){
					int a=possible_hits[i];double b=possible_distance[i];
					for(int j=0;j<i;j++){
						if(possible_distance[j]>b){
							possible_distance.erase(possible_distance.begin()+i);
							possible_hits.erase(possible_hits.begin()+i);
							possible_distance.insert(possible_distance.begin()+j,b);
							possible_hits.insert(possible_hits.begin()+j,a);
						}
					}
				}
				
				//check through subsequent detectors
				for(int i=0;(unsigned)i<possible_hits.size();i++){
					if(valid_dets[part][possible_hits[i]]){//if this detector is valid
						for(int j=0;j<=i;j++){//count through all before
							if(i==j){
								det_hits[part].push_back(possible_hits[i]);
								stop_distance[part]=possible_distance[i];
							}else if(!valid_doubles[possible_hits[i]][possible_hits[j]]){//is this before one paired
								i=possible_hits.size();//finish if not
								j=i;//finish if not
							}
						}
					}else{
						if(i==0)stop_distance[part]=possible_distance[i];//save if no hits
						i=possible_hits.size();//finish if it hits an impassable
					}
				}
			}//obstructed_sub
		}else{stop_distance[part]=0;}//momentum zero if
	}//paticle validity if
}	




//check all detectors and all particles
void exp_core::det_check_hits_all(bool obstruct){
	/*if(!decay_events)*/this->det_check_hits(0,obstruct);
	if(!targ_fusion){//save on pointless calculation
		this->det_check_hits(1,obstruct);		
	}	
	
	
	//decay events more complex	
	if(decay_events){
		TVector3 tA(0,0,0),tB(0,0,0);
		
		//if there is a lifetime event which may need implant adjustment
		if(use_offsets&&lifetime>0&&implant_mode!=0){
			
			if(!obstruct)this->max_distance_hit(0);
			//calculate the implant position, if on non-obstruct mode it hasn't been calculated yet.
			
			if(stop_distance[0]>0){//if its -1 no hit/no implant OR if 0 never left so offset decay will be correct
				if(stop_distance[0]<(offset_decay-offset_primary).Mag()){//if implanted before decay
					offset_decay-=offset_primary;
					offset_decay.SetMag(stop_distance[0]);//change the decay position
					offset_decay+=offset_primary;
					
					if(post_beta>0){//reverse the boost to the lab frame for decay particles
						TVector3 boost=-lorR.Vect();
						if(boost.Mag()>0){
							boost.SetMag(post_beta);
							lorA.Boost(boost);
							lorB.Boost(boost);
						}
					}
				
					tA=lorA.Vect();//create offset vectors for implant hit calcs
					if(tA.Mag()>0){tA*=implant_mode;tA.SetMag(0.01);}
					
					tB=lorB.Vect();
					if(tB.Mag()>0){tB*=implant_mode;tB.SetMag(0.01);}
				}
			}
		}
		
		offset_decay+=(tA);//shift JUST before or after implant detector (if needed)
		this->det_check_hits(2,obstruct);
		offset_decay-=(tA);//return it to normal
		
		offset_decay+=(tB);
		if(decay_type==1) this->det_check_hits(3,false);//gamma never obstructed
		else this->det_check_hits(3,obstruct);
		offset_decay-=(tB);
	}
}
	
	
	
	
int exp_core::current_multiplicity(){
	int mult=0;
	for(int i=0;i<4;i++){//loop over 4 particles count multiplicity
		if(det_hits[i].size()>0)mult++;
	}
	return mult;
}
	



//check all detectors and all particles
double exp_core::max_distance_hit(int part){
	if(stop_distance[part]<0){
		double length=0;	
		
		for(int j=0;(unsigned)j<det_hits[part].size();j++){
			double zxy=detectors[det_hits[part][j]].hit_pos_vec(lor_point[part],offsets[part]).Mag();
			if(length<zxy){length=zxy;}
		}

		if(length>0){
			stop_distance[part]=length;
			return length;
		}
		return -1;
	}
	return stop_distance[part];
}

/*

//call AFTER det_check_hits_all
void exp_core::implant_correct(bool obstruct){
	
	if(decay_events&&use_offsets){
		
		this->max_distance_hit(0);
		if(stop_distance[0]>=0){
		
			
			if(stop_distance[0]<(offset_decay-offset_primary).Mag()){
				offset_decay-=offset_primary;
				offset_decay.SetMag(stop_distance[0]);
				offset_decay+=offset_primary;
				
				if(post_beta>0){
					TVector3 boost=-lorR.Vect();
					if(boost.Mag()>0){
						boost.SetMag(post_beta);
						lorA.Boost(boost);
						lorB.Boost(boost);
					}
				}		
			}
			
			
			TVector3 temp=lorA.Vect();
			if(temp.Mag()>0){
				temp.SetMag(0.01);offset_decay-=temp;
				this->det_check_hits(2,obstruct);
				offset_decay+=temp;
			}
			
			temp=lorB.Vect();
			if(temp.Mag()>0){
				temp.SetMag(0.01);offset_decay-=temp;
				
				if(decay_type==1) this->det_check_hits(3,false);//gamma never obstructed
				else this->det_check_hits(3,obstruct);
				
				offset_decay+=temp;
			}
		}
	}
}*/


//add all hits to their detectors run after check hits
void exp_core::addhitsdet(){
	/*if(!decay_events)*/for(int i=0;(unsigned)i<det_hits[0].size();i++){
		if(use_offsets)detectors[det_hits[0][i]].Fill(lor_point[0],offsets[0]);else detectors[det_hits[0][i]].Fill(lor_point[0]);
	}
	if(!targ_fusion){//save on pointless calculation
		for(int i=0;(unsigned)i<det_hits[1].size();i++){if(use_offsets)detectors[det_hits[1][i]].Fill(lor_point[1],offsets[1]);else detectors[det_hits[1][i]].Fill(lor_point[1]);}		
	}	
	if(decay_events){
		for(int i=0;(unsigned)i<det_hits[2].size();i++){if(use_offsets)detectors[det_hits[2][i]].Fill(lor_point[2],offsets[2]);else detectors[det_hits[2][i]].Fill(lor_point[2]);}
		for(int i=0;(unsigned)i<det_hits[3].size();i++){if(use_offsets)detectors[det_hits[3][i]].Fill(lor_point[3],offsets[3]);else detectors[det_hits[3][i]].Fill(lor_point[3]);}
	}
}

void exp_core::extend_doubles()
{
	//set rows
	for(int i=valid_doubles.size();(unsigned)i<detectors.size();i++)valid_doubles.push_back(vector<bool>(valid_doubles.size(),false));
	
	//extend column
	for(int i=0;(unsigned)i<valid_doubles.size();i++)
		for(int j=valid_doubles[i].size();(unsigned)j<detectors.size();j++)if(i==j)valid_doubles[i].push_back(true);else valid_doubles[i].push_back(false);
}