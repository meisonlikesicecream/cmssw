#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletConfigBuilder_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletConfigBuilder_h

#include <vector>
#include <utility>
#include <set>
#include <iostream>
#include <fstream>
#include <cstdlib>

namespace trklet {

  class TrackletConfigBuilder{

  public:
  
    TrackletConfigBuilder(bool combinedmodules);
    
    void writeAll(std::ostream& wires, std::ostream& memories, std::ostream& modules);
  
  private:
    
    std::pair<unsigned int,unsigned int> seedLayers(unsigned int iSeed);
    
    void initGeom();
    
    std::pair<double,double> seedRadii(unsigned int iseed);
    
    bool validTEPair(unsigned int iseed, unsigned int iTE1, unsigned int iTE2);
    
    void buildTE();
    
    void buildTC();
    
    std::pair<double, double> seedPhiRange(double rproj, unsigned int iSeed, unsigned int iTC);
    
    void buildProjections();
    
    double phi(double r1, double phi1, double r2, double phi2,double r);
    
    double rinv(double r1, double phi1, double r2, double phi2);
    
    std::string iSeedStr(unsigned int iSeed);
    
    std::string numStr(unsigned int i);
    
    std::string iTCStr(unsigned int iTC);
    
    std::string iRegStr(unsigned int iReg, unsigned int iSeed);
    
    std::string TCName(unsigned int iSeed,unsigned int iTC);
    
    std::string LayerName(unsigned int ilayer);
    
    std::string TPROJName(unsigned int iSeed, unsigned int iTC,
			  unsigned int ilayer, unsigned int ireg);
    
    std::string PRName(unsigned int ilayer, unsigned int ireg);
    
    void writeProjectionMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    std::string SPName(unsigned int l1, unsigned int ireg1, unsigned int ivm1,
		       unsigned int l2, unsigned int ireg2, unsigned int ivm2,
		       unsigned int iseed);
    
    std::string TEName(unsigned int l1, unsigned int ireg1, unsigned int ivm1,
		       unsigned int l2, unsigned int ireg2, unsigned int ivm2,
		       unsigned int iseed);
    
    std::string TCNAme(unsigned int iseed, unsigned int iTC);
    
    void writeSPMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    void writeAPMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    void writeCMMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    void writeVMPROJMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    void writeFMMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    void writeASMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    void writeVMSMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);

    void writeTPARMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);

    void writeTFMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    void writeCTMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    void writeILMemories(std::ostream& os, std::ostream& memories, std::ostream& modules);
    
    unsigned int NSector_; //Number of sectors
    double rcrit_; //critical radius that defines the sector
    
    bool combinedmodules_; //if true write configuration for combined modules
    
    double rinvmax_; //Max value for valid rinv
    double rmaxdisk_; //Maximim disk radius
    double zlength_; //Maximim (abslute) z-positon in barrel
    double rmean_[6];  //Mean layer radius  
    double zmean_[5];  //Mean disk z-position  
    
    double dphisectorHG_; //Full sector width
    
    unsigned int NTC_[8]; //Number of TC per seeding combination

    unsigned int NTPSeedRegion_[8]; //Number of TP per region (outer) 
  
    unsigned int NRegions_[11]; //Regions (all stubs memories 6 layers +5 disks
    unsigned int NVMME_[11]; //Number of MEs (all stubs memories 6 layers +5 disks
    std::pair<unsigned int, unsigned int> NVMTE_[8]; //number of TEs for each seeding combination
    

    std::vector<std::pair<double,double> > allStubs_[11]; //Allstubs phi ranges
    std::vector<std::pair<double,double> > VMStubsME_[11]; //ME VMstubs phi ranges
    std::pair<std::vector<std::pair<double, double > >,
      std::vector<std::pair<double, double > > > VMStubsTE_[8]; //TE VMstubs phi ranges
    
    std::vector<std::pair<unsigned int ,unsigned int> > TE_[8];
    
    std::vector<std::vector<unsigned int> > TC_[8];
    
    std::vector< std::vector<std::pair<unsigned int, unsigned int> > > projections_[11]; //seedindex and TC
    
  };
}
#endif
