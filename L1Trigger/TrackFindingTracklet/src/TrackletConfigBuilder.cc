#include <vector>
#include <utility>
#include <set>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>

#include "L1Trigger/TrackFindingTracklet/interface/TrackletConfigBuilder.h"

using namespace std;
using namespace trklet;

TrackletConfigBuilder::TrackletConfigBuilder(bool combinedmodules){
  
  NSector_=9;
  rcrit_=55.0;

  combinedmodules_=combinedmodules;
  
  rinvmax_=0.01*0.3*3.8/2.0; //0.01 to convert to cm-1
  
  rmaxdisk_=120.0;
  zlength_=120.0;
    
  rmean_[0]=(rmaxdisk_*851)/4096;
  rmean_[1]=(rmaxdisk_*1269)/4096;
  rmean_[2]=(rmaxdisk_*1784)/4096;
  rmean_[3]=(rmaxdisk_*2347)/4096;
  rmean_[4]=(rmaxdisk_*2936)/4096;
  rmean_[5]=(rmaxdisk_*3697)/4096;
  
  zmean_[0]=(zlength_*2239)/2048;
  zmean_[1]=(zlength_*2645)/2048;
  zmean_[2]=(zlength_*3163)/2048;
  zmean_[3]=(zlength_*3782)/2048;
  zmean_[4]=(zlength_*4523)/2048;

  double rsectmin=21.8;
  double rsectmax=112.7;
  
  dphisectorHG_ = 2 * M_PI / NSector_ +
      rinvmax_*std::max(rcrit_ - rsectmin, rsectmax - rcrit_);
  
  NRegions_[0]=8;
  NRegions_[1]=4;
  NRegions_[2]=4;
  NRegions_[3]=4;
  NRegions_[4]=4;
  NRegions_[5]=4;
  NRegions_[6]=4;
  NRegions_[7]=4;
  NRegions_[8]=4;
  NRegions_[9]=4;
  NRegions_[10]=4;
  
  NVMME_[0]=4;
  NVMME_[1]=8;
  NVMME_[2]=8;
  NVMME_[3]=8;
  NVMME_[4]=8;
  NVMME_[5]=8;
  NVMME_[6]=8;
  NVMME_[7]=4;
  NVMME_[8]=4;
  NVMME_[9]=4;
  NVMME_[10]=4;

  
  NTPSeedRegion_[0]=1;
  NTPSeedRegion_[1]=1;
  NTPSeedRegion_[2]=1;
  NTPSeedRegion_[3]=1;
  NTPSeedRegion_[4]=1;
  NTPSeedRegion_[5]=1;
  NTPSeedRegion_[6]=1;
  NTPSeedRegion_[7]=1;

  NVMTE_[0]=std::pair<unsigned int, unsigned int>(4,8);
  NVMTE_[1]=std::pair<unsigned int, unsigned int>(4,4);
  NVMTE_[2]=std::pair<unsigned int, unsigned int>(4,8);
  NVMTE_[3]=std::pair<unsigned int, unsigned int>(4,8);
  NVMTE_[4]=std::pair<unsigned int, unsigned int>(4,4);
  NVMTE_[5]=std::pair<unsigned int, unsigned int>(4,4);
  NVMTE_[6]=std::pair<unsigned int, unsigned int>(2,4);
  NVMTE_[7]=std::pair<unsigned int, unsigned int>(2,4);
  
  initGeom();
    
  buildTE();
  
  NTC_[0]=12;
  NTC_[1]=4;
  NTC_[2]=4;  //8
  NTC_[3]=4;
  NTC_[4]=4;
  NTC_[5]=4;
  NTC_[6]=8;  //4
  NTC_[7]=4;
    
  buildTC();

  buildProjections();

}

std::pair<unsigned int,unsigned int> TrackletConfigBuilder::seedLayers(unsigned int iSeed) {

  if (iSeed==0) return std::pair<unsigned int, unsigned int>(0,1);
  if (iSeed==1) return std::pair<unsigned int, unsigned int>(1,2);
  if (iSeed==2) return std::pair<unsigned int, unsigned int>(2,3);
  if (iSeed==3) return std::pair<unsigned int, unsigned int>(4,5);
  if (iSeed==4) return std::pair<unsigned int, unsigned int>(6,7);
  if (iSeed==5) return std::pair<unsigned int, unsigned int>(8,9);
  if (iSeed==6) return std::pair<unsigned int, unsigned int>(0,6);
  
  assert(iSeed==7);
  
  return std::pair<unsigned int, unsigned int>(1,6);
  
}


void TrackletConfigBuilder::initGeom(){
  
  for(unsigned int ilayer=0;ilayer<11;ilayer++) {
    double dphi= dphisectorHG_/NRegions_[ilayer];
    for(unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++){
      std::vector< std::pair<unsigned int, unsigned int> > emptyVec;
      projections_[ilayer].push_back(emptyVec);
      double phimin=dphi*iReg;
      double phimax=phimin+dphi;
      std::pair<double,double> tmp(phimin,phimax);
      allStubs_[ilayer].push_back(tmp);
      double dphiVM=dphi/NVMME_[ilayer];
      for(unsigned int iVM=0;iVM<NVMME_[ilayer];iVM++){
	double phivmmin=phimin+iVM*dphiVM;
	double phivmmax=phivmmin+dphiVM;
	std::pair<double,double> tmp(phivmmin,phivmmax);
	VMStubsME_[ilayer].push_back(tmp);
      }
    }
  }
  for(unsigned int iseed=0;iseed<8;iseed++){
    unsigned int l1=seedLayers(iseed).first;
    unsigned int l2=seedLayers(iseed).second;
    unsigned int nVM1=NVMTE_[iseed].first;
    unsigned int nVM2=NVMTE_[iseed].second;
    double dphiVM=dphisectorHG_/(nVM1*NRegions_[l1]);
    for(unsigned int iVM=0;iVM<nVM1*NRegions_[l1];iVM++){
      double phivmmin=iVM*dphiVM;
      double phivmmax=phivmmin+dphiVM;
      std::pair<double,double> tmp(phivmmin,phivmmax);
      VMStubsTE_[iseed].first.push_back(tmp);
    }
    dphiVM=dphisectorHG_/(nVM2*NRegions_[l2]);
    for(unsigned int iVM=0;iVM<nVM2*NRegions_[l2];iVM++){
      double phivmmin=iVM*dphiVM;
      double phivmmax=phivmmin+dphiVM;
      std::pair<double,double> tmp(phivmmin,phivmmax);
      VMStubsTE_[iseed].second.push_back(tmp);
    }
  }    
}

std::pair<double,double> TrackletConfigBuilder::seedRadii(unsigned int iseed) {
  
  std::pair<unsigned int,unsigned int> seedlayers=seedLayers(iseed);
  
  unsigned int l1=seedlayers.first;
  unsigned int l2=seedlayers.second;
  
  double r1,r2;
  
  if (iseed<4) {  //barrel seeding
    r1=rmean_[l1];
    r2=rmean_[l2];
  } else if (iseed<6) { //disk seeding
    r1=rmean_[0]+40.0; //Somwwhat of a hack - but allows finding all the regions
    //when projecting to L1
    r2=r1*zmean_[l2-6]/zmean_[l1-6];
  } else { //overlap seeding
    r1=rmean_[l1];
    r2=r1*zmean_[l2-6]/zlength_;
  }
  
  return std::pair<double,double>(r1,r2);
}


bool TrackletConfigBuilder::validTEPair(unsigned int iseed, unsigned int iTE1, unsigned int iTE2){
  
  double rinvmin=999.9;
  double rinvmax=-999.9;
  
  double phi1[2]={VMStubsTE_[iseed].first[iTE1].first,VMStubsTE_[iseed].first[iTE1].second};
  double phi2[2]={VMStubsTE_[iseed].second[iTE2].first,VMStubsTE_[iseed].second[iTE2].second};
  
  std::pair<double, double> seedradii=seedRadii(iseed);
  
  for (unsigned int i1=0;i1<2;i1++) {
    for (unsigned int i2=0;i2<2;i2++) {
      double arinv=rinv(seedradii.first,phi1[i1],seedradii.second,phi2[i2]);
      if (arinv<rinvmin) rinvmin=arinv;
      if (arinv>rinvmax) rinvmax=arinv;
    }
  }
  
  if (rinvmin>rinvmax_) return false;
  if (rinvmax<-rinvmax_) return false;
  
  return true;
  
}
  
void TrackletConfigBuilder::buildTE() {
  
  for(unsigned int iseed=0;iseed<8;iseed++){
    for(unsigned int i1=0;i1<VMStubsTE_[iseed].first.size();i1++){
      for(unsigned int i2=0;i2<VMStubsTE_[iseed].second.size();i2++){
	if (validTEPair(iseed,i1,i2)) {
	  std::pair<unsigned int, unsigned int> tmp(i1,i2);
	  TE_[iseed].push_back(tmp);
	}
      }
    }
  } 
}


void TrackletConfigBuilder::buildTC(){

  for(unsigned int iSeed=0;iSeed<8;iSeed++){
    unsigned int nTC=NTC_[iSeed];
    std::vector<std::pair<unsigned int, unsigned int> >& TEs=TE_[iSeed];
    std::vector<std::vector<unsigned int> >& TCs=TC_[iSeed];
    
    //Very naive method to group TEs in TC
    
    double invnTC=nTC*(1.0/TEs.size());

    for(unsigned int iTE=0;iTE<TEs.size();iTE++){
      int iTC=invnTC*iTE;
      assert(iTC<(int)nTC);
      if (iTC>=(int)TCs.size()) {
	std::vector<unsigned int> tmp;
	tmp.push_back(iTE);
	TCs.push_back(tmp);
      } else {
	TCs[iTC].push_back(iTE);
      }
    }
  }
}

std::pair<double, double> TrackletConfigBuilder::seedPhiRange(double rproj, unsigned int iSeed, unsigned int iTC){
  
  std::vector<std::vector<unsigned int> >& TCs=TC_[iSeed];
  
  std::pair<double,double> seedradii=seedRadii(iSeed);
  
  double phimin=999.0;
  double phimax=-999.0;
  for (unsigned int iTE=0;iTE<TCs[iTC].size();iTE++) {
    unsigned int theTE=TCs[iTC][iTE];
    unsigned int l1TE=TE_[iSeed][theTE].first;
    unsigned int l2TE=TE_[iSeed][theTE].second;
    double phi1[2]={VMStubsTE_[iSeed].first[l1TE].first,
		    VMStubsTE_[iSeed].first[l1TE].second};
    double phi2[2]={VMStubsTE_[iSeed].second[l2TE].first,
		    VMStubsTE_[iSeed].second[l2TE].second};
    for(unsigned int i1=0;i1<2;i1++){
      for(unsigned int i2=0;i2<2;i2++){
	double aphi=phi(seedradii.first, phi1[i1],seedradii.second, phi2[i2],
			rproj);
	if (aphi<phimin) phimin=aphi;
	if (aphi>phimax) phimax=aphi;
      }
    }
  }
  return std::pair<double, double>(phimin,phimax);
}

void TrackletConfigBuilder::buildProjections(){
  
  //FIXME - should be member data
  
  //                       L1 L2 L3 L4 L5 L6 D1 D2 D3 D4 D5 
  int matchport[8][11]= { {-1,-1, 1, 2, 3, 4, 4, 3, 2, 1,-1},   //L1L2
			  { 1,-1,-1, 2, 3,-1, 4, 3, 2, 1,-1},   //L2L3
			  { 1, 2,-1,-1, 3, 4, 4, 3,-1,-1,-1},   //L3L4
			  { 1, 2, 3, 4,-1,-1,-1,-1,-1,-1,-1},   //L5L6
			  { 1, 2,-1,-1,-1,-1,-1,-1, 2, 3, 4},   //D1D2
			  { 1,-1,-1,-1,-1,-1, 2, 3,-1,-1, 4},   //D3D4
			  {-1,-1,-1,-1,-1,-1,-1, 1, 2, 3, 4},   //L1D1
			  { 1,-1,-1,-1,-1,-1,-1, 2, 3, 4,-1}};  //L2D1
  
  for(unsigned int iseed=0;iseed<8;iseed++){ 
    std::vector<std::vector<unsigned int> >& TCs=TC_[iseed];
    
    for (unsigned int ilayer=0;ilayer<11;ilayer++) {
      if (matchport[iseed][ilayer]==-1) continue;
      for (unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++) {
	for (unsigned int iTC=0;iTC<TCs.size();iTC++){
	  double rproj=rmaxdisk_;
	  if (ilayer<6) rproj=rmean_[ilayer];
	  std::pair<double, double> phiRange=seedPhiRange(rproj,iseed,iTC);
	  if (phiRange.first<allStubs_[ilayer][iReg].second&&
	      phiRange.second>allStubs_[ilayer][iReg].first) {
	    std::pair<unsigned int, unsigned int> tmp(iseed,iTC); //seedindex and TC
	    projections_[ilayer][iReg].push_back(tmp);
	  }
	}
      }
    }
  }    
}

  
double TrackletConfigBuilder::phi(double r1, double phi1, double r2, double phi2,double r) {
  
  double rhoinv=rinv(r1,phi1,r2,phi2);
  if (fabs(rhoinv)>rinvmax_) {
    rhoinv=rinvmax_*rhoinv/fabs(rhoinv);
    }
  return phi1+asin(0.5*r*rhoinv)-asin(0.5*r1*rhoinv);
  
}


    
double TrackletConfigBuilder::rinv(double r1, double phi1, double r2, double phi2) {
  
  double deltaphi=phi1-phi2;
  return 2*sin(deltaphi)/sqrt(r2*r2+r1*r1-2*r1*r2*cos(deltaphi));
  
}

std::string TrackletConfigBuilder::iSeedStr(unsigned int iSeed) {

  static std::string name[8]={"L1L2","L2L3","L3L4","L5L6","D1D2","D3D4","L1D1","L2D1"};
  
  assert(iSeed<8);
  return name[iSeed];
  
}


std::string TrackletConfigBuilder::numStr(unsigned int i) {
  
  static std::string num[32]={"1","2","3","4","5","6","7","8","9","10",
			      "11","12","13","14","15","16","17","18","19","20",
			      "21","22","23","24","25","26","27","28","29","30","31","32"};
  assert(i<32);
  return num[i];
  
}

std::string TrackletConfigBuilder::iTCStr(unsigned int iTC) {
  
  static std::string name[12]={"A","B","C","D","E","F","G","H","I","J","K","L"};
  
  assert(iTC<12);
  return name[iTC];
  
}

std::string TrackletConfigBuilder::iRegStr(unsigned int iReg, unsigned int iSeed) {
  
  static std::string name[8]={"A","B","C","D","E","F","G","H"};
  
  static std::string nameOverlap[8]={"X","Y","Z","W","Q","R","S","T"};
  
  static std::string nameL2L3[4]={"I","J","K","L"};
  
  if (iSeed==1) {
    assert(iReg<4);
    return nameL2L3[iReg];
  }
  if (iSeed==6||iSeed==7) {
    assert(iReg<8);
    return nameOverlap[iReg];
  }
  assert(iReg<8);
  return name[iReg];
  
}



std::string TrackletConfigBuilder::TCName(unsigned int iSeed,unsigned int iTC){
  if (combinedmodules_) {
    return "TP_"+iSeedStr(iSeed)+iTCStr(iTC);
  } else {
    return "TC_"+iSeedStr(iSeed)+iTCStr(iTC);
  }
}

std::string TrackletConfigBuilder::LayerName(unsigned int ilayer) {
  return ilayer<6?("L"+numStr(ilayer)):("D"+numStr(ilayer-6));
}

std::string TrackletConfigBuilder::TPROJName(unsigned int iSeed, unsigned int iTC,
					     unsigned int ilayer, unsigned int ireg) {
  return "TPROJ_"+iSeedStr(iSeed)+iTCStr(iTC)+"_"+LayerName(ilayer)+"PHI"+iTCStr(ireg);
}

std::string TrackletConfigBuilder::PRName(unsigned int ilayer, unsigned int ireg) {
  if (combinedmodules_) {
    return "MP_"+LayerName(ilayer)+"PHI"+iTCStr(ireg);
  } else {
    return "PR_"+LayerName(ilayer)+"PHI"+iTCStr(ireg);
  }
}

void TrackletConfigBuilder::writeProjectionMemories(std::ostream& os, std::ostream& memories, std::ostream&){
  
  for(unsigned int ilayer=0;ilayer<11;ilayer++){
    for(unsigned int ireg=0;ireg<projections_[ilayer].size();ireg++){
      for(unsigned int imem=0;imem<projections_[ilayer][ireg].size();imem++){
	
	unsigned int iSeed=projections_[ilayer][ireg][imem].first;
	unsigned int iTC=projections_[ilayer][ireg][imem].second;
	
	memories << "TrackletProjections: "+TPROJName(iSeed,iTC,ilayer,ireg)+" [54]"<<std::endl;
	
	os<<TPROJName(iSeed,iTC,ilayer,ireg)
	  <<" input=> "<<TCName(iSeed,iTC)<<".projout"<<LayerName(ilayer)<<"PHI"<<iTCStr(ireg)
	  <<" output=> "<<PRName(ilayer,ireg)<<".projin"<<std::endl;
	
      }	
    }
  }    
}


std::string TrackletConfigBuilder::SPName(unsigned int l1, unsigned int ireg1, unsigned int ivm1,
					  unsigned int l2, unsigned int ireg2, unsigned int ivm2,
					  unsigned int iseed) {
  
  return "SP_"+LayerName(l1)+"PHI"+iRegStr(ireg1,iseed)+numStr(ivm1)+
    "_"+LayerName(l2)+"PHI"+iRegStr(ireg2,iseed)+numStr(ivm2);
    
}

std::string TrackletConfigBuilder::TEName(unsigned int l1, unsigned int ireg1, unsigned int ivm1,
					  unsigned int l2, unsigned int ireg2, unsigned int ivm2,
					  unsigned int iseed) {
  
  return "TE_"+LayerName(l1)+"PHI"+iRegStr(ireg1,iseed)+numStr(ivm1)+
    "_"+LayerName(l2)+"PHI"+iRegStr(ireg2,iseed)+numStr(ivm2);
  
}

std::string TrackletConfigBuilder::TCNAme(unsigned int iseed, unsigned int iTC){
  return "TC_"+iSeedStr(iseed)+iTCStr(iTC);
}

void TrackletConfigBuilder::writeSPMemories(std::ostream& os, std::ostream& memories, std::ostream& modules){

  if (combinedmodules_) return;
  
  for(unsigned int iSeed=0;iSeed<8;iSeed++) {
    
    for(unsigned int iTC=0;iTC<TC_[iSeed].size();iTC++){
      for(unsigned int iTE=0;iTE<TC_[iSeed][iTC].size();iTE++){
	
	unsigned int theTE=TC_[iSeed][iTC][iTE];
	
	unsigned int TE1=TE_[iSeed][theTE].first;
	unsigned int TE2=TE_[iSeed][theTE].second;

	unsigned int l1=seedLayers(iSeed).first;
	unsigned int l2=seedLayers(iSeed).second;
	
	memories << "StubPairs: "
		 << SPName(l1,TE1/NVMTE_[iSeed].first,TE1,l2,TE2/NVMTE_[iSeed].second,TE2,iSeed)
		 << " [12]"<<std::endl;
	modules << "TrackletEngine: "
		<< TEName(l1,TE1/NVMTE_[iSeed].first,TE1,l2,TE2/NVMTE_[iSeed].second,TE2,iSeed)
		<< std::endl;
	
	os << SPName(l1,TE1/NVMTE_[iSeed].first,TE1,l2,TE2/NVMTE_[iSeed].second,TE2,iSeed)
	   <<	" input=> "<<TEName(l1,TE1/NVMTE_[iSeed].first,TE1,l2,TE2/NVMTE_[iSeed].second,TE2,iSeed)
	   << ".stubpairout output=> "<<TCNAme(iSeed,iTC)<<".stubpairin"<<std::endl;
      }
    }
  }
}


void TrackletConfigBuilder::writeAPMemories(std::ostream& os, std::ostream& memories, std::ostream& modules){

  if (combinedmodules_) return;
  
  for(unsigned int ilayer=0;ilayer<11;ilayer++){
    for(unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++){
      
      memories << "AllProj: AP_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<" [56]"<<std::endl;
      modules << "ProjectionRouter: PR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<std::endl;
      
      os << "AP_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	 << " input=> PR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	 << ".allprojout output=> MC_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	 << ".allprojin"<<std::endl;
      
    }
  }
}

void TrackletConfigBuilder::writeCMMemories(std::ostream& os, std::ostream& memories, std::ostream& modules){

  if (combinedmodules_) return;
  
  for(unsigned int ilayer=0;ilayer<11;ilayer++){
    for(unsigned int iME=0;iME<NVMME_[ilayer]*NRegions_[ilayer];iME++){
      
      memories << "CandidateMatch: CM_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iME/NVMME_[ilayer])
	       <<iME+1<<" [12]"<<std::endl;
      modules << "MatchEngine: ME_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iME/NVMME_[ilayer])
	      <<iME+1<<std::endl;
      
      os << "CM_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iME/NVMME_[ilayer])<<iME+1
	 << " input=> ME_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iME/NVMME_[ilayer])<<iME+1
	 << ".matchout output=> MC_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iME/NVMME_[ilayer])
	 << ".matchin"<<std::endl;
      
    }
  }
  
}


void TrackletConfigBuilder::writeVMPROJMemories(std::ostream& os, std::ostream& memories, std::ostream&){

  if (combinedmodules_) return;
  
  for(unsigned int ilayer=0;ilayer<11;ilayer++){      
    for(unsigned int iME=0;iME<NVMME_[ilayer]*NRegions_[ilayer];iME++){
      
      memories << "VMProjections: VMPROJ_"<<LayerName(ilayer)<<"PHI"
	       <<iTCStr(iME/NVMME_[ilayer])<<iME+1<<" [13]"<<std::endl;
      
      os << "VMPROJ_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iME/NVMME_[ilayer])<<iME+1
	 << " input=> PR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iME/NVMME_[ilayer])
	 << ".vmprojout"<<"PHI"<<iTCStr(iME/NVMME_[ilayer])<<iME+1<<" output=> ME_"
	 << LayerName(ilayer)<<"PHI"<<iTCStr(iME/NVMME_[ilayer])
	 << iME+1 << ".vmprojin"<<std::endl;
    }
  }

}
 

void TrackletConfigBuilder::writeFMMemories(std::ostream& os, std::ostream& memories, std::ostream& modules){
  //                       L1 L2 L3 L4 L5 L6 D1 D2 D3 D4 D5 
  int matchport[8][11]= { {-1,-1, 1, 2, 3, 4, 4, 3, 2, 1,-1},   //L1L2
			  { 1,-1,-1, 2, 3,-1, 4, 3, 2, 1,-1},   //L2L3
			  { 1, 2,-1,-1, 3, 4, 4, 3,-1,-1,-1},   //L3L4
			  { 1, 2, 3, 4,-1,-1,-1,-1,-1,-1,-1},   //L5L6
			  { 1, 2,-1,-1,-1,-1,-1,-1, 2, 3, 4},   //D1D2
			  { 1,-1,-1,-1,-1,-1, 2, 3,-1,-1, 4},   //D3D4
			  {-1,-1,-1,-1,-1,-1,-1, 1, 2, 3, 4},   //L1D1
			  { 1,-1,-1,-1,-1,-1,-1, 2, 3, 4,-1}};  //L2D1


  if (combinedmodules_) {
    for(unsigned int ilayer=0;ilayer<11;ilayer++){
      for(unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++){
	modules << "MatchProcessor: MP_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<std::endl;
	for(unsigned int iSeed=0;iSeed<8;iSeed++) {
	  if (matchport[iSeed][ilayer]==-1) continue;
	  memories << "FullMatch: FM_"<<iSeedStr(iSeed)<<"_"<<LayerName(ilayer)<<"PHI"
		   <<iTCStr(iReg)<<" [36]"<<std::endl;
	  os << "FM_"<<iSeedStr(iSeed)<<"_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	     << " input=> MP_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	     << ".matchout1 output=> FT_"<<iSeedStr(iSeed)<<".fullmatch"
	     << matchport[iSeed][ilayer]<<"in"<<iReg+1<<std::endl;
	}
      }
    }	
  } else {
    for(unsigned int ilayer=0;ilayer<11;ilayer++){
      for(unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++){
	modules << "MatchCalculator: MC_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<std::endl;
	for(unsigned int iSeed=0;iSeed<8;iSeed++) {
	  if (matchport[iSeed][ilayer]==-1) continue;
	  memories << "FullMatch: FM_"<<iSeedStr(iSeed)<<"_"<<LayerName(ilayer)<<"PHI"
		   <<iTCStr(iReg)<<" [36]"<<std::endl;
	  os << "FM_"<<iSeedStr(iSeed)<<"_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	     << " input=> MC_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	     << ".matchout1 output=> FT_"<<iSeedStr(iSeed)<<".fullmatch"
	     << matchport[iSeed][ilayer]<<"in"<<iReg+1<<std::endl;
	}
      }
    }	
  }      
}

void TrackletConfigBuilder::writeASMemories(std::ostream& os, std::ostream& memories, std::ostream& modules){

  if (combinedmodules_) {
   //First write AS memories used by MatchProcessor
    for(unsigned int ilayer=0;ilayer<11;ilayer++){
      for(unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++){
	memories << "AllStubs: AS_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<"n1"<<" [42]"<<std::endl;
	if (combinedmodules_) {
	  modules << "VMRouterCM: VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<std::endl;
	} else {
	  modules << "VMRouter: VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<std::endl;
	}
	os << "AS_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<"n1"
	   << " input=> VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	   << ".allstubout output=> MP_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	   << ".allstubin"<<std::endl;
      }
    }
    
    
    //Next write AS memories used by TrackletProcessor
    for(unsigned int ilayer=0;ilayer<11;ilayer++){
      for(int iReg=0;iReg<(int)NRegions_[ilayer];iReg++){
	
	unsigned int nmem=1;
	
	for(unsigned int iSeed=0;iSeed<8;iSeed++) {
	  
	  unsigned int l1=seedLayers(iSeed).first;
	  unsigned int l2=seedLayers(iSeed).second;
	  
	  if (ilayer!=l1&&ilayer!=l2) continue;
	  
	  for(unsigned int iTC=0;iTC<TC_[iSeed].size();iTC++){

	    int nTCReg=TC_[iSeed].size()/NRegions_[l2];
	    
	    int iTCReg=iTC/nTCReg;
            
	    int jTCReg=iTC%nTCReg;
	    
	    if (ilayer==l2) {
	      if (iTCReg!=iReg)
		continue;
	    }

	    string ext="";

	    if (ilayer==l1) {
	      int ratio=NRegions_[l1]/NRegions_[l2];
	      //std::cout << "ratio : "<<ratio << std::endl;
	      int min=iTCReg*ratio-1+jTCReg;
	      int max=(iTCReg+1)*ratio-(nTCReg-jTCReg-1);
	      if ((int)iReg<min || (int)iReg>max)
		continue;

	      if (max-min>=2) {
		ext="_M";
		if (iReg==min)
		  ext="_R";
		if (iReg==max)
		  ext="_L";
	      }
	      
	      if (max-min==1) {
		if (nTCReg==2) {
		  if (jTCReg==0){
		    if (iReg==min)
		      ext="_R";
		    if (iReg==max)
		      ext="_B";
		  }
		  if (jTCReg==1){
		    if (iReg==min)
		      ext="_A";
		    if (iReg==max)
		      ext="_L";
		  }
		}
		if (nTCReg==3) {
		  if (jTCReg==0){
		    if (iReg==min)
		      ext="_R";
		    if (iReg==max)
		      ext="_F";
		  }
		  if (jTCReg==1){
		    if (iReg==min)
		      ext="_E";
		    if (iReg==max)
		      ext="_D";
		  }		  
		  if (jTCReg==2){
		    if (iReg==min)
		      ext="_C";
		    if (iReg==max)
		      ext="_L";
		  }		  
		}
	      }
	      //cout << "iseed ratio max min : "<<iSeed<<" "<<ratio<<" "<<max<<" "<<min<<endl;
	      assert(ext!="");
	      
	    }

	    nmem++;
	    memories << "AllStubs: AS_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<"n"<<nmem
		     <<" [42]"<<std::endl;
	    os << "AS_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<"n"<<nmem
	       << " input=> VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	       << ".allstubout"<<ext<<" output=> TP_"<<iSeedStr(iSeed)<<iTCStr(iTC);
	    if (ilayer==l1){
	      os << ".innerallstubin"<<std::endl;
	    } else {
	      os << ".outerallstubin"<<std::endl;
	    }
	      
	  }
	}
      }
    }
  } else {
    //First write AS memories used by MatchCalculator
    for(unsigned int ilayer=0;ilayer<11;ilayer++){
      for(unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++){
	memories << "AllStubs: AS_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<"n1"<<" [42]"<<std::endl;
	if (combinedmodules_) {
	  modules << "VMRouterCM: VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<std::endl;
	} else {
	  modules << "VMRouter: VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<std::endl;
	}
	os << "AS_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<"n1"
	   << " input=> VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	   << ".allstubout output=> MC_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	   << ".allstubin"<<std::endl;
      }
    }
    
    
    //Next write AS memories used by TrackletCalculator
    for(unsigned int ilayer=0;ilayer<11;ilayer++){
      for(unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++){
	
	unsigned int nmem=1;
	
	for(unsigned int iSeed=0;iSeed<8;iSeed++) {
	  
	  unsigned int l1=seedLayers(iSeed).first;
	  unsigned int l2=seedLayers(iSeed).second;
	  
	  if (ilayer!=l1&&ilayer!=l2) continue;
	  
	  for(unsigned int iTC=0;iTC<TC_[iSeed].size();iTC++){
	    
	    bool used=false;
	    for(unsigned int iTE=0;iTE<TC_[iSeed][iTC].size();iTE++){
	      
	      unsigned int theTE=TC_[iSeed][iTC][iTE];
	      
	      unsigned int TE1=TE_[iSeed][theTE].first;
	      unsigned int TE2=TE_[iSeed][theTE].second;
	      
	      if (l1==ilayer && iReg==TE1/NVMTE_[iSeed].first) used=true;
	      if (l2==ilayer && iReg==TE2/NVMTE_[iSeed].second) used=true;
	    }
	    
	    if (used) {
	      nmem++;
	      memories << "AllStubs: AS_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<"n"<<nmem
		       <<" [42]"<<std::endl;
	      os << "AS_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<"n"<<nmem
		 << " input=> VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
		 << ".allstubout output=> TC_"<<iSeedStr(iSeed)<<iTCStr(iTC);
	      if (ilayer==l1){
		os << ".innerallstubin"<<std::endl;
	      } else {
		os << ".outerallstubin"<<std::endl;
	      }
	      
	    }
	  }
	}
      }
    }
  }
}

void TrackletConfigBuilder::writeVMSMemories(std::ostream& os, std::ostream& memories, std::ostream&){

  if (combinedmodules_) {
    //First write VMS memories used by MatchProcessor
    for(unsigned int ilayer=0;ilayer<11;ilayer++){
      for(unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++){
	memories << "VMStubsME: VMSME_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
		 <<"n1 [18]"<<std::endl;
	os << "VMSME_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<"n1"
	   << " input=> VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	   << ".vmstuboutPHI"<<iTCStr(iReg)
	   <<" output=> MP_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	   << ".vmstubin"<<std::endl;
      }
    }
  
    //Next write VMS memories used by TrackletProcessor
    for(unsigned int iSeed=0;iSeed<8;iSeed++) {
      
      //FIXME - code could be cleaner
      unsigned int l1=seedLayers(iSeed).first;
      unsigned int l2=seedLayers(iSeed).second;
      
      unsigned int ilayer=seedLayers(iSeed).second;
      
      //for(unsigned int iReg=0;iReg<NRegions_[ilayer];iReg++){
      
      unsigned int nTCReg=TC_[iSeed].size()/NRegions_[l2];

      for(unsigned int iReg=0;iReg<NRegions_[l2];iReg++){

	unsigned int nmem=0;
	//Hack since we use same module twice
	if (iSeed==7) nmem=2;

	for (unsigned iTC=0;iTC<nTCReg;iTC++) {
	
	  nmem++;
	  memories << "VMStubsTE: VMSTE_"<<LayerName(ilayer)<<"PHI"<<iRegStr(iReg,iSeed)
		   <<"n"<<nmem<<" [18]"<<std::endl;
	  os << "VMSTE_"<<LayerName(ilayer)<<"PHI"<<iRegStr(iReg,iSeed)<<"n"<<nmem
	     << " input=> VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	     << ".vmstubout_seed_"<<iSeed
	     << " output=> TP_"<<LayerName(l1)<<LayerName(l2)<<iTCStr(iReg*nTCReg+iTC)
	     << ".outervmstubin"<<std::endl;
	}
      }
    }
  } else {
    //First write VMS memories used by MatchEngine
    for(unsigned int ilayer=0;ilayer<11;ilayer++){
      for(unsigned int iVMME=0;iVMME<NVMME_[ilayer]*NRegions_[ilayer];iVMME++){
	unsigned int iReg=iVMME/NVMME_[ilayer];
	memories << "VMStubsME: VMSME_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<iVMME+1
		 <<"n1 [18]"<<std::endl;
	os << "VMSME_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<iVMME+1<<"n1"
	   << " input=> VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	   << ".vmstuboutPHI"<<iTCStr(iReg)<<iVMME+1
	   <<" output=> ME_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)<<iVMME+1
	   << ".vmstubin"<<std::endl;
      }
    }
  
    //Next write VMS memories used by MatchCalculator
    for(unsigned int iSeed=0;iSeed<8;iSeed++) {
      
      for (unsigned int innerouterseed=0;innerouterseed<2;innerouterseed++){
	
	//FIXME - code could be cleaner
	unsigned int l1=seedLayers(iSeed).first;
	unsigned int l2=seedLayers(iSeed).second;
      
	unsigned int NVMTE1=NVMTE_[iSeed].first;
	unsigned int NVMTE2=NVMTE_[iSeed].second;
	
      
	unsigned int ilayer=seedLayers(iSeed).first;
	unsigned int NVMTE=NVMTE_[iSeed].first;
	if (innerouterseed==1) {
	  ilayer=seedLayers(iSeed).second;
	  NVMTE=NVMTE_[iSeed].second;
	}
      
	for(unsigned int iVMTE=0;iVMTE<NVMTE*NRegions_[ilayer];iVMTE++){
	  unsigned int iReg=iVMTE/NVMTE;
	
	  unsigned int nmem=0;
	
	  for(unsigned int iTE=0;iTE<TE_[iSeed].size();iTE++){

	    unsigned int TE1=TE_[iSeed][iTE].first;
	    unsigned int TE2=TE_[iSeed][iTE].second;
	  
	    bool used=false;
	  
	    if ( innerouterseed==0 && iVMTE==TE1) used=true;
	    if ( innerouterseed==1 && iVMTE==TE2) used=true;
	  
	    if (!used) continue;
	    nmem++;
	    memories << "VMStubsTE: VMSTE_"<<LayerName(ilayer)<<"PHI"<<iRegStr(iReg,iSeed)<<iVMTE+1
		     <<"n"<<nmem<<" [18]"<<std::endl;
	    os << "VMSTE_"<<LayerName(ilayer)<<"PHI"<<iRegStr(iReg,iSeed)<<iVMTE+1<<"n"<<nmem
	       << " input=> VMR_"<<LayerName(ilayer)<<"PHI"<<iTCStr(iReg)
	       << ".vmstuboutPHI"<<iTCStr(iReg)<<iVMTE+1
	       << " output=> TE_"<<LayerName(l1)<<"PHI"<<iRegStr(TE1/NVMTE1,iSeed)<<TE1+1
	       << "_"<<LayerName(l2)<<"PHI"<<iRegStr(TE2/NVMTE2,iSeed)<<TE2+1;
	    if (innerouterseed==0) {
	      os << ".innervmstubin"<<std::endl;
	    } else {
	      os << ".outervmstubin"<<std::endl;
	    }
	  }
	}
      }
    }
  }
}

void TrackletConfigBuilder::writeTPARMemories(std::ostream& os, std::ostream& memories, std::ostream& modules){

  if (combinedmodules_) {
    for(unsigned int iSeed=0;iSeed<8;iSeed++) {
      for(unsigned int iTP=0;iTP<TC_[iSeed].size();iTP++){
	memories << "TrackletParameters: TPAR_"<<iSeedStr(iSeed)<<iTCStr(iTP)
	       <<" [56]"<<std::endl;
	modules << "TrackletProcessor: TP_"<<iSeedStr(iSeed)<<iTCStr(iTP)<<std::endl;
	os << "TPAR_"<<iSeedStr(iSeed)<<iTCStr(iTP)
	   << " input=> TP_"<<iSeedStr(iSeed)<<iTCStr(iTP)
	   << ".trackpar output=> FT_"<<iSeedStr(iSeed)<<".tparin"<<std::endl;
      }
    }
  } else {
    for(unsigned int iSeed=0;iSeed<8;iSeed++) {
      for(unsigned int iTC=0;iTC<TC_[iSeed].size();iTC++){
	memories << "TrackletParameters: TPAR_"<<iSeedStr(iSeed)<<iTCStr(iTC)
	       <<" [56]"<<std::endl;
	modules << "TrackletCalculator: TC_"<<iSeedStr(iSeed)<<iTCStr(iTC)<<std::endl;
	os << "TPAR_"<<iSeedStr(iSeed)<<iTCStr(iTC)
	   << " input=> TC_"<<iSeedStr(iSeed)<<iTCStr(iTC)
	   << ".trackpar output=> FT_"<<iSeedStr(iSeed)<<".tparin"<<std::endl;
      }
    }
  }
}

void TrackletConfigBuilder::writeTFMemories(std::ostream& os, std::ostream& memories, std::ostream& modules){
  
  for(unsigned int iSeed=0;iSeed<8;iSeed++) {
    memories << "TrackFit: TF_"<<iSeedStr(iSeed)<<" [126]"<<std::endl;
    modules << "FitTrack: FT_"<<iSeedStr(iSeed)<<std::endl;
    os << "TF_"<<iSeedStr(iSeed)
       << " input=> FT_"<<iSeedStr(iSeed)
       << ".trackout output=> PD.trackin"<<std::endl;
  }
  
}

void TrackletConfigBuilder::writeCTMemories(std::ostream& os, std::ostream& memories, std::ostream& modules){
  
  modules << "PurgeDuplicate: PD"<<std::endl;
  
    for(unsigned int iSeed=0;iSeed<8;iSeed++) {
      memories << "CleanTrack: CT_"<<iSeedStr(iSeed)<<" [126]"<<std::endl;
      os << "CT_"<<iSeedStr(iSeed)
	 << " input=> PD.trackout output=>"<<std::endl;
    }
    
}

void TrackletConfigBuilder::writeILMemories(std::ostream& os, std::ostream& memories, std::ostream&){

  //FIXME these should not be hardcoded - but for now wanted to avoid reading file
  string dtcname[52];
  unsigned int layerdisk[52];
  double phimin[52];
  double phimax[52];

  //FIXME layerdisk numbering should be 0 to 10
  
  dtcname[0] ="PS10G_1";    layerdisk[0] =1;   phimin[0] = 0.304273;   phimax[0] = 0.742925;
  dtcname[1] ="PS10G_1";    layerdisk[1] =7;   phimin[1] = -0.185672;  phimax[1] =  0.883803;
  dtcname[2] ="PS10G_1";    layerdisk[2] =9;   phimin[2] = -0.132414;  phimax[2] =  0.830545;
  dtcname[3] ="PS10G_1";    layerdisk[3] =11;  phimin[3] = -0.132414;  phimax[3] =  0.830545;
  dtcname[4] ="PS10G_2";    layerdisk[4] =1;   phimin[4] = -0.0133719; phimax[4] =  0.715599;
  dtcname[5] ="PS10G_2";    layerdisk[5] =8;   phimin[5] = -0.110089;  phimax[5] =  0.808221;
  dtcname[6] ="PS10G_2";    layerdisk[6] =10;  phimin[6] = -0.132414;  phimax[6] =  0.830545;
  dtcname[7] ="PS10G_3";    layerdisk[7] = 2;  phimin[7] = -0.11381;   phimax[7] =  0.822812;
  dtcname[8] ="PS10G_3";    layerdisk[8] = 8;  phimin[8] = -0.185672;  phimax[8] =  0.883803;
  dtcname[9] ="PS10G_4";    layerdisk[9] = 7;  phimin[9] = -0.0823971; phimax[9] =  0.780529;
  dtcname[10]="PS10G_4";    layerdisk[10]= 9;  phimin[10]= -0.0963091; phimax[10]=  0.794441;
  dtcname[11]="PS10G_4";    layerdisk[11]= 11; phimin[11]= -0.0963091; phimax[11]=  0.794441;
  dtcname[12]="PS_1";       layerdisk[12]= 3;  phimin[12]= 0.0827748;  phimax[12]=  0.615357;
  dtcname[13]="PS_1";       layerdisk[13]= 8;  phimin[13]= -0.0823971; phimax[13]=  0.780529;
  dtcname[14]="PS_2";       layerdisk[14]= 3;  phimin[14]= -0.0917521; phimax[14]=  0.614191;
  dtcname[15]="PS_2";       layerdisk[15]= 10; phimin[15]= -0.0963091; phimax[15]=  0.794441;
  dtcname[16]="2S_1";       layerdisk[16]= 4;  phimin[16]= -0.0246209; phimax[16]=  0.763311;
  dtcname[17]="2S_1";       layerdisk[17]= 5;  phimin[17]= 0.261875;   phimax[17]=  0.403311;
  dtcname[18]="2S_2";       layerdisk[18]= 5;  phimin[18]= -0.0542445; phimax[18]=  0.715509;
  dtcname[19]="2S_3";       layerdisk[19]= 6;  phimin[19]= 0.0410126;  phimax[19]=  0.730605;
  dtcname[20]="2S_4";       layerdisk[20]= 6;  phimin[20]= -0.0428961; phimax[20]=  0.693862;
  dtcname[21]="2S_4";       layerdisk[21]= 9;  phimin[21]= -0.0676705; phimax[21]=  0.765802;
  dtcname[22]="2S_5";       layerdisk[22]= 7;  phimin[22]= -0.0648206; phimax[22]=  0.762952;
  dtcname[23]="2S_5";       layerdisk[23]= 10; phimin[23]= -0.0676705; phimax[23]=  0.765802;
  dtcname[24]="2S_6";       layerdisk[24]= 8;  phimin[24]= -0.0648206; phimax[24]=  0.762952;
  dtcname[25]="2S_6";       layerdisk[25]= 11; phimin[25]= -0.0676705; phimax[25]=  0.765802;
  dtcname[26]="negPS10G_1"; layerdisk[26]= 1;  phimin[26]= -0.023281;  phimax[26]=  0.372347;
  dtcname[27]="negPS10G_1"; layerdisk[27]= 7;  phimin[27]= -0.185672;  phimax[27]=  0.883803;
  dtcname[28]="negPS10G_1"; layerdisk[28]= 9;  phimin[28]= -0.132414;  phimax[28]=  0.830545;
  dtcname[29]="negPS10G_1"; layerdisk[29]= 11; phimin[29]= -0.132414;  phimax[29]=  0.830545;
  dtcname[30]="negPS10G_2"; layerdisk[30]= 1;  phimin[30]= -0.0133719; phimax[30]=  0.715599;
  dtcname[31]="negPS10G_2"; layerdisk[31]= 8;  phimin[31]= -0.110089;  phimax[31]=  0.808221;
  dtcname[32]="negPS10G_2"; layerdisk[32]= 10; phimin[32]= -0.132414;  phimax[32]=  0.830545;
  dtcname[33]="negPS10G_3"; layerdisk[33]= 2;  phimin[33]= -0.115834;  phimax[33]=  0.813823;
  dtcname[34]="negPS10G_3"; layerdisk[34]= 8;  phimin[34]= -0.185672;  phimax[34]=  0.883803;
  dtcname[35]="negPS10G_4"; layerdisk[35]= 7;  phimin[35]= -0.0823971; phimax[35]=  0.780529;
  dtcname[36]="negPS10G_4"; layerdisk[36]= 9;  phimin[36]= -0.0963091; phimax[36]=  0.794441;
  dtcname[37]="negPS10G_4"; layerdisk[37]= 11; phimin[37]= -0.0963091; phimax[37]=  0.794441;
  dtcname[38]="negPS_1";    layerdisk[38]= 3;  phimin[38]= -0.0961318; phimax[38]=  0.445198;
  dtcname[39]="negPS_1";    layerdisk[39]= 8;  phimin[39]= -0.0823971; phimax[39]=  0.780529;
  dtcname[40]="negPS_2";    layerdisk[40]= 3;  phimin[40]= -0.0917521; phimax[40]=  0.614191;
  dtcname[41]="negPS_2";    layerdisk[41]= 10; phimin[41]= -0.0963091; phimax[41]=  0.794441;
  dtcname[42]="neg2S_1";    layerdisk[42]= 4;  phimin[42]= -0.0246209; phimax[42]=  0.763311;
  dtcname[43]="neg2S_1";    layerdisk[43]= 5;  phimin[43]= 0.261875;   phimax[43]=  0.403311;
  dtcname[44]="neg2S_2";    layerdisk[44]= 5;  phimin[44]= -0.0542445; phimax[44]=  0.715509;
  dtcname[45]="neg2S_3";    layerdisk[45]= 6;  phimin[45]= 0.0410126;  phimax[45]=  0.730605;
  dtcname[46]="neg2S_4";    layerdisk[46]= 6;  phimin[46]= -0.0428961; phimax[46]=  0.693862;
  dtcname[47]="neg2S_4";    layerdisk[47]= 9;  phimin[47]= -0.06767;   phimax[47]=  0.765802;
  dtcname[48]="neg2S_5";    layerdisk[48]= 7;  phimin[48]= -0.0648201; phimax[48]=  0.762952;
  dtcname[49]="neg2S_5";    layerdisk[49]= 10; phimin[49]= -0.06767;   phimax[49]=  0.765802;
  dtcname[50]="neg2S_6";    layerdisk[50]= 8;  phimin[50]= -0.0648201; phimax[50]=  0.762952;
  dtcname[51]="neg2S_6";    layerdisk[51]= 11; phimin[51]= -0.06767;   phimax[51]=  0.765802;
  
   
  double dphi=0.5*dphisectorHG_-M_PI/NSector_;
  
  for(unsigned int i=0;i<52;i++) {
    
    double phimintmp=phimin[i]+dphi;
    double phimaxtmp=phimax[i]+dphi;
      
    for(unsigned int iReg=0; iReg<NRegions_[layerdisk[i]-1]; iReg++) {
      
      if (allStubs_[layerdisk[i]-1][iReg].first>phimaxtmp &&
	  allStubs_[layerdisk[i]-1][iReg].second<phimintmp) continue;

      if (allStubs_[layerdisk[i]-1][iReg].second<phimaxtmp) { 
	memories << "InputLink: IL_"<<LayerName(layerdisk[i]-1)<<"PHI"
		 <<iTCStr(iReg)<<"_"<<dtcname[i]<<"_A"<< " [36]"<<std::endl;
	os << "IL_"<<LayerName(layerdisk[i]-1)<<"PHI"<<iTCStr(iReg)<<"_"<<dtcname[i]<<"_A"
	   << " input=> output=> VMR_"<<LayerName(layerdisk[i]-1)<<"PHI"<<iTCStr(iReg)
	   << ".stubin"<<std::endl;
      }
      
      if (allStubs_[layerdisk[i]-1][iReg].first>phimintmp) { 
	memories << "InputLink: IL_"<<LayerName(layerdisk[i]-1)<<"PHI"
		 <<iTCStr(iReg)<<"_"<<dtcname[i]<<"_B"<< " [36]"<<std::endl;
	os << "IL_"<<LayerName(layerdisk[i]-1)<<"PHI"<<iTCStr(iReg)<<"_"<<dtcname[i]<<"_B"
	   << " input=> output=> VMR_"<<LayerName(layerdisk[i]-1)<<"PHI"<<iTCStr(iReg)
	   << ".stubin"<<std::endl;
      }
      
    }
  }
  
}

void TrackletConfigBuilder::writeAll(std::ostream& wires, std::ostream& memories, std::ostream& modules){

  writeILMemories(wires,memories,modules);
  writeASMemories(wires,memories,modules);
  writeVMSMemories(wires,memories,modules);
  writeSPMemories(wires,memories,modules);
  writeProjectionMemories(wires,memories,modules);
  writeTPARMemories(wires,memories,modules);
  writeVMPROJMemories(wires,memories,modules);
  writeAPMemories(wires,memories,modules);
  writeCMMemories(wires,memories,modules);
  writeFMMemories(wires,memories,modules);
  writeTFMemories(wires,memories,modules);
  writeCTMemories(wires,memories,modules);

}



