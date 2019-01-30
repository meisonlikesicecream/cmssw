// File: SiStripDigitizerAlgorithm.cc
// Description:  Steering class for digitization.
// Modified 15/May/2013 mark.grimes@bristol.ac.uk - Modified so that the digi-sim link has the correct
// index for the sim hits stored. It was previously always set to zero (I won't mention that it was
// me who originally wrote that).
// Modified on Feb 11, 2015: prolay.kumar.mal@cern.ch & Jean-Laurent.Agram@cern.ch
//                           Added/Modified the individual strip noise in zero suppression
//                           mode from the conditions DB; previously, the digitizer used to
//                           consider the noise value for individual strips inside a module from
//                           the central strip noise value.
//////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <algorithm>
#include <iostream>
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SiStripDigitizerAlgorithm.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/GeometrySurface/interface/BoundSurface.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "CalibTracker/Records/interface/SiStripDependentRecords.h"
#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
#include "CondFormats/SiStripObjects/interface/SiStripThreshold.h"
#include "CondFormats/SiStripObjects/interface/SiStripPedestals.h"
#include "CondFormats/SiStripObjects/interface/SiStripBadStrip.h"
#include "CLHEP/Random/RandFlat.h"

#include <string.h>
#include <sstream>

#include <boost/algorithm/string.hpp>

SiStripDigitizerAlgorithm::SiStripDigitizerAlgorithm(const edm::ParameterSet& conf):
  lorentzAngleName(conf.getParameter<std::string>("LorentzAngle")),
  theThreshold(conf.getParameter<double>("NoiseSigmaThreshold")),
  cmnRMStib(conf.getParameter<double>("cmnRMStib")),
  cmnRMStob(conf.getParameter<double>("cmnRMStob")),
  cmnRMStid(conf.getParameter<double>("cmnRMStid")),
  cmnRMStec(conf.getParameter<double>("cmnRMStec")),
  APVSaturationProbScaling_(conf.getParameter<double>("APVSaturationProbScaling")),
  makeDigiSimLinks_(conf.getUntrackedParameter<bool>("makeDigiSimLinks", false)),
  peakMode(conf.getParameter<bool>("APVpeakmode")),
  noise(conf.getParameter<bool>("Noise")),
  RealPedestals(conf.getParameter<bool>("RealPedestals")), 
  SingleStripNoise(conf.getParameter<bool>("SingleStripNoise")),
  CommonModeNoise(conf.getParameter<bool>("CommonModeNoise")),
  BaselineShift(conf.getParameter<bool>("BaselineShift")),
  APVSaturationFromHIP(conf.getParameter<bool>("APVSaturationFromHIP")),
  theFedAlgo(conf.getParameter<int>("FedAlgorithm")),
  zeroSuppression(conf.getParameter<bool>("ZeroSuppression")),
  theElectronPerADC(conf.getParameter<double>( peakMode ? "electronPerAdcPeak" : "electronPerAdcDec" )),

  apv_minAmplitude(conf.getParameter<double>( "apv_minAmplitude" )), // 10000;
  apv_decayConstantInMicroS(conf.getParameter<double>( "apv_decayConstantInMicroS" )), // 200;
  apv_nPreviousInteractionsToSimulate(conf.getParameter<unsigned int>( "apv_nPreviousInteractionsToSimulate" )), //  = 100;
  apv_nBaselineToGenerate(conf.getParameter<unsigned int>( "apv_nBaselineToGenerate" )), // = 2000;
  apv_smallChangeThreshold(conf.getParameter<double>( "apv_smallChangeThreshold" )), // 0.001 
  apv_smallChangeN(conf.getParameter<unsigned int>( "apv_smallChangeN" )), //5

  apv_maxResponse(conf.getParameter<double>( "apv_maxResponse" )), // 729

  apv_rate(conf.getParameter<double>( "apv_rate" )), // 66.2
  apv_mVPerQ(conf.getParameter<double>( "apv_mVPerQ" )), // 5.5
  apv_fCPerElectron(conf.getParameter<double>( "apvfCPerElectron" )), // 1.602e-4

  theTOFCutForPeak(conf.getParameter<double>("TOFCutForPeak")),
  theTOFCutForDeconvolution(conf.getParameter<double>("TOFCutForDeconvolution")),
  tofCut(peakMode ? theTOFCutForPeak : theTOFCutForDeconvolution),
  cosmicShift(conf.getUntrackedParameter<double>("CosmicDelayShift")),
  inefficiency(conf.getParameter<double>("Inefficiency")),
  pedOffset((unsigned int)conf.getParameter<double>("PedestalsOffset")),
  PreMixing_(conf.getParameter<bool>("PreMixingMode")),
  theSiHitDigitizer(new SiHitDigitizer(conf)),
  theSiPileUpSignals(new SiPileUpSignals()),
  theSiNoiseAdder(new SiGaussianTailNoiseAdder(theThreshold)),
  theSiDigitalConverter(new SiTrivialDigitalConverter(theElectronPerADC, PreMixing_)),
  theSiZeroSuppress(new SiStripFedZeroSuppression(theFedAlgo)),
  APVProbabilityFile(conf.getParameter<edm::FileInPath>("APVProbabilityFile")) {

  if (peakMode) {
    LogDebug("StripDigiInfo")<<"APVs running in peak mode (poor time resolution)";
  } else {
    LogDebug("StripDigiInfo")<<"APVs running in deconvolution mode (good time resolution)";
  };
  if(SingleStripNoise) LogDebug("SiStripDigitizerAlgorithm")<<" SingleStripNoise: ON";
  else LogDebug("SiStripDigitizerAlgorithm")<<" SingleStripNoise: OFF";
  if(CommonModeNoise) LogDebug("SiStripDigitizerAlgorithm")<<" CommonModeNoise: ON";
  else LogDebug("SiStripDigitizerAlgorithm")<<" CommonModeNoise: OFF";
  if(APVSaturationFromHIP){  
    std::string line; 
    APVProbaFile.open((APVProbabilityFile.fullPath()).c_str());
    if (APVProbaFile.is_open()){
      while ( getline (APVProbaFile,line) ){
        std::vector<std::string> strs;
          boost::split(strs,line,boost::is_any_of(" "));
          if(strs.size()==2){
            mapOfAPVprobabilities[std::stoi(strs.at(0))]=std::stof(strs.at(1));
          }
      }
      APVProbaFile.close();
    }else throw cms::Exception("MissingInput")
         << "It seems that the APV probability list is missing\n";
  }

  apvBaselineDistributions_tib_.push_back( TH1F("TIB1", "TIB1", 100, 0, 729) );
  apvBaselineDistributions_tib_.push_back( TH1F("TIB2", "TIB2", 100, 0, 729) );
  apvBaselineDistributions_tib_.push_back( TH1F("TIB3", "TIB3", 100, 0, 729) );
  apvBaselineDistributions_tib_.push_back( TH1F("TIB4", "TIB4", 100, 0, 729) );


  apvBaselineDistributions_tob_.push_back( TH1F("TOB1", "TOB1", 100, 0, 729) );
  apvBaselineDistributions_tob_.push_back( TH1F("TOB2", "TOB2", 100, 0, 729) );
  apvBaselineDistributions_tob_.push_back( TH1F("TOB3", "TOB3", 100, 0, 729) );
  apvBaselineDistributions_tob_.push_back( TH1F("TOB4", "TOB4", 100, 0, 729) );
  apvBaselineDistributions_tob_.push_back( TH1F("TOB5", "TOB5", 100, 0, 729) );
  apvBaselineDistributions_tob_.push_back( TH1F("TOB6", "TOB6", 100, 0, 729) );


  apvBaselineDistributions_tid_.push_back( TH1F("TID1", "TID1", 100, 0, 729) );
  apvBaselineDistributions_tid_.push_back( TH1F("TID2", "TID2", 100, 0, 729) );
  apvBaselineDistributions_tid_.push_back( TH1F("TID3", "TID3", 100, 0, 729) );

  apvBaselineDistributions_tec_.push_back( TH1F("TEC1", "TEC1", 100, 0, 729) );
  apvBaselineDistributions_tec_.push_back( TH1F("TEC2", "TEC2", 100, 0, 729) );
  apvBaselineDistributions_tec_.push_back( TH1F("TEC3", "TEC3", 100, 0, 729) );
  apvBaselineDistributions_tec_.push_back( TH1F("TEC4", "TEC4", 100, 0, 729) );
  apvBaselineDistributions_tec_.push_back( TH1F("TEC5", "TEC5", 100, 0, 729) );
  apvBaselineDistributions_tec_.push_back( TH1F("TEC6", "TEC6", 100, 0, 729) );
  apvBaselineDistributions_tec_.push_back( TH1F("TEC7", "TEC7", 100, 0, 729) );
}

SiStripDigitizerAlgorithm::~SiStripDigitizerAlgorithm(){
}

void
SiStripDigitizerAlgorithm::initializeDetUnit(StripGeomDetUnit const * det, const edm::EventSetup& iSetup){
  edm::ESHandle<SiStripBadStrip> deadChannelHandle;
  iSetup.get<SiStripBadChannelRcd>().get(deadChannelHandle);

  unsigned int detId = det->geographicalId().rawId();
  int numStrips = (det->specificTopology()).nstrips();  

  SiStripBadStrip::Range detBadStripRange = deadChannelHandle->getRange(detId);
  //storing the bad strip of the the module. the module is not removed but just signal put to 0
  std::vector<bool>& badChannels = allBadChannels[detId];
  badChannels.clear();
  badChannels.insert(badChannels.begin(), numStrips, false);
  for(SiStripBadStrip::ContainerIterator it = detBadStripRange.first; it != detBadStripRange.second; ++it) {
    SiStripBadStrip::data fs = deadChannelHandle->decode(*it);
    for(int strip = fs.firstStrip; strip < fs.firstStrip + fs.range; ++strip) {
		badChannels[strip] = true;
	}
  }
  firstChannelsWithSignal[detId] = numStrips;
  lastChannelsWithSignal[detId]= 0;

  //  if(APVSaturationFromHIP){
  //  std::bitset<6> &bs=SiStripTrackerAffectedAPVMap[detId];
  //  if(bs.any())theAffectedAPVvector.push_back(std::make_pair(detId,bs));
  //}

}

void
SiStripDigitizerAlgorithm::initializeEvent(const edm::EventSetup& iSetup) {
  theSiPileUpSignals->reset();
  // This should be clear by after all calls to digitize(), but I might as well make sure
  associationInfoForDetId_.clear();

  APVSaturationProb_ = APVSaturationProbScaling_;  // reset probability
  SiStripTrackerAffectedAPVMap.clear();
  FirstLumiCalc_ = true;
  FirstDigitize_ = true;

  //get gain noise pedestal lorentzAngle from ES handle
  edm::ESHandle<ParticleDataTable> pdt;
  iSetup.getData(pdt);
  setParticleDataTable(&*pdt);
  iSetup.get<SiStripLorentzAngleSimRcd>().get(lorentzAngleName,lorentzAngleHandle);

  for ( unsigned int i = 0; i < apvBaselineDistributions_tib_.size(); ++ i ) apvBaselineDistributions_tib_[i].Reset();
  for ( unsigned int i = 0; i < apvBaselineDistributions_tob_.size(); ++ i ) apvBaselineDistributions_tob_[i].Reset();
  for ( unsigned int i = 0; i < apvBaselineDistributions_tid_.size(); ++ i ) apvBaselineDistributions_tid_[i].Reset();
  for ( unsigned int i = 0; i < apvBaselineDistributions_tec_.size(); ++ i ) apvBaselineDistributions_tec_[i].Reset();

}

//  Run the algorithm for a given module
//  ------------------------------------

void
SiStripDigitizerAlgorithm::accumulateSimHits(std::vector<PSimHit>::const_iterator inputBegin,
                                             std::vector<PSimHit>::const_iterator inputEnd,
                                             size_t inputBeginGlobalIndex,
					     unsigned int tofBin,
                                             const StripGeomDetUnit* det,
                                             const GlobalVector& bfield,
					     const TrackerTopology *tTopo,
                                             CLHEP::HepRandomEngine* engine) {
  // produce SignalPoints for all SimHits in detector
  unsigned int detID = det->geographicalId().rawId();
  int numStrips = (det->specificTopology()).nstrips();  

  size_t thisFirstChannelWithSignal = numStrips;
  size_t thisLastChannelWithSignal = 0;

  float langle = (lorentzAngleHandle.isValid()) ? lorentzAngleHandle->getLorentzAngle(detID) : 0.;

  std::vector<float> locAmpl(numStrips, 0.);

  // Loop over hits

  uint32_t detId = det->geographicalId().rawId();
  // First: loop on the SimHits
  if(CLHEP::RandFlat::shoot(engine) > inefficiency) {
    AssociationInfoForChannel* pDetIDAssociationInfo; // I only need this if makeDigiSimLinks_ is true...
    if( makeDigiSimLinks_ ) pDetIDAssociationInfo=&(associationInfoForDetId_[detId]); // ...so only search the map if that is the case
    std::vector<float> previousLocalAmplitude; // Only used if makeDigiSimLinks_ is true. Needed to work out the change in amplitude.

    size_t simHitGlobalIndex=inputBeginGlobalIndex; // This needs to stored to create the digi-sim link later
    for (std::vector<PSimHit>::const_iterator simHitIter = inputBegin; simHitIter != inputEnd; ++simHitIter, ++simHitGlobalIndex ) {
      // skip hits not in this detector.
      if((*simHitIter).detUnitId() != detId) {
        continue;
      }
      // check TOF
      if (std::fabs(simHitIter->tof() - cosmicShift - det->surface().toGlobal(simHitIter->localPosition()).mag()/30.) < tofCut && simHitIter->energyLoss()>0) {
        if( makeDigiSimLinks_ ) previousLocalAmplitude=locAmpl; // Not needed except to make the sim link association.
        size_t localFirstChannel = numStrips;
        size_t localLastChannel  = 0;
        // process the hit
        theSiHitDigitizer->processHit(&*simHitIter, *det, bfield, langle, locAmpl, localFirstChannel, localLastChannel, tTopo, engine);
 
        if(thisFirstChannelWithSignal > localFirstChannel) thisFirstChannelWithSignal = localFirstChannel;
        if(thisLastChannelWithSignal < localLastChannel) thisLastChannelWithSignal = localLastChannel;

        if( makeDigiSimLinks_ ) { // No need to do any of this if truth association was turned off in the configuration
          for( size_t stripIndex=0; stripIndex<locAmpl.size(); ++stripIndex ) {
            // Work out the amplitude from this SimHit from the difference of what it was before and what it is now
            float signalFromThisSimHit=locAmpl[stripIndex]-previousLocalAmplitude[stripIndex];
            if( signalFromThisSimHit!=0 ) { // If this SimHit had any contribution I need to record it.
              auto& associationVector=(*pDetIDAssociationInfo)[stripIndex];
              bool addNewEntry=true;
              // Make sure the hit isn't in already. I've seen this a few times, it always seems to happen in pairs so I think
              // it's something to do with the stereo strips.
              for( auto& associationInfo : associationVector ) {
                if( associationInfo.trackID==simHitIter->trackId() && associationInfo.eventID==simHitIter->eventId() ) {
                  // The hit is already in, so add this second contribution and move on
                  associationInfo.contributionToADC+=signalFromThisSimHit;
                  addNewEntry=false;
                  break;
                }
              } // end of loop over associationVector
              // If the hit wasn't already in create a new association info structure.
              if( addNewEntry ) associationVector.push_back( AssociationInfo{ simHitIter->trackId(), simHitIter->eventId(), signalFromThisSimHit, simHitGlobalIndex, tofBin } );
            } // end of "if( signalFromThisSimHit!=0 )"
          } // end of loop over locAmpl strips
        } // end of "if( makeDigiSimLinks_ )"
      } // end of TOF check
    } // end for
  }
  theSiPileUpSignals->add(detID, locAmpl, thisFirstChannelWithSignal, thisLastChannelWithSignal);

  if(firstChannelsWithSignal[detID] > thisFirstChannelWithSignal) firstChannelsWithSignal[detID] = thisFirstChannelWithSignal;
  if(lastChannelsWithSignal[detID] < thisLastChannelWithSignal) lastChannelsWithSignal[detID] = thisLastChannelWithSignal;
}

//============================================================================                
void SiStripDigitizerAlgorithm::calculateInstlumiScale(PileupMixingContent* puInfo){
  //Instlumi scalefactor calculating for dynamic inefficiency                                 

  if (puInfo && FirstLumiCalc_) {

    const std::vector<int> bunchCrossing = puInfo->getMix_bunchCrossing();
    const std::vector<float> TrueInteractionList = puInfo->getMix_TrueInteractions();
    const int bunchSpacing = puInfo->getMix_bunchSpacing();                                 

    double RevFreq = 11245.;
    double minBXsec = 70.0E-27;  // use 70mb as an approximation
    double Bunch = 2100.;        // 2016 value
    if (bunchSpacing == 50) Bunch = Bunch/2.;

    int pui = 0, p = 0;
    std::vector<int>::const_iterator pu;
    std::vector<int>::const_iterator pu0 = bunchCrossing.end();

    for (pu=bunchCrossing.begin(); pu!=bunchCrossing.end(); ++pu) {
      if (*pu==0) {
        pu0 = pu;
        p = pui;
      }
      pui++;
    }
    if (pu0!=bunchCrossing.end()) {  // found the in-time interaction
      double Tintr = TrueInteractionList.at(p);
      double instLumi = Bunch*Tintr*RevFreq/minBXsec;
      APVSaturationProb_ = instLumi/6.0E33;      
    }
    FirstLumiCalc_ = false;
  }
}

//============================================================================                

void SiStripDigitizerAlgorithm::calcuateAPVBaselines(
              TrackingGeometry::DetUnitContainer detUnits,
              const TrackerTopology *tTopo
              ) {

    std::vector<unsigned int> nStrips_tib = { 0, 0, 0, 0 };
    std::vector<unsigned int> nStrips_tob = { 0, 0, 0, 0, 0, 0 };
    std::vector<unsigned int> nStrips_tid = { 0, 0, 0 };
    std::vector<unsigned int> nStrips_tec = { 0, 0, 0, 0, 0, 0, 0 };

    std::vector<unsigned int> nStripsNonZero_tib = { 0, 0, 0, 0 };
    std::vector<unsigned int> nStripsNonZero_tob = { 0, 0, 0, 0, 0, 0 };
    std::vector<unsigned int> nStripsNonZero_tid = { 0, 0, 0 };
    std::vector<unsigned int> nStripsNonZero_tec = { 0, 0, 0, 0, 0, 0, 0 };

    std::vector<float> occupancy_tib = { 0, 0, 0, 0 };
    std::vector<float> occupancy_tob = { 0, 0, 0, 0, 0, 0 };
    std::vector<float> occupancy_tid = { 0, 0, 0 };
    std::vector<float> occupancy_tec = { 0, 0, 0, 0, 0, 0, 0 };

    TH1F h_charge_tib = TH1F("h_charge_tib", "h_charge_tib", 60, 0, 60000);
    TH1F h_charge_tob = TH1F("h_charge_tob", "h_charge_tob", 60, 0, 60000);
    TH1F h_charge_tid = TH1F("h_charge_tid", "h_charge_tid", 60, 0, 60000);
    TH1F h_charge_tec = TH1F("h_charge_tec", "h_charge_tec", 60, 0, 60000);

    for(TrackingGeometry::DetUnitContainer::const_iterator iu = detUnits.begin(); iu != detUnits.end(); iu ++){

      auto sgd = dynamic_cast<StripGeomDetUnit const*>((*iu));
      if (sgd != 0){

        DetId detId( sgd->geographicalId() );
        uint32_t SubDet = detId.subdetId();
        int layer = -1;

        if(SubDet==3) {
          layer = tTopo->tibLayer(detId);
          nStrips_tib[layer-1] += sgd->specificTopology().nstrips();
        }
        else if(SubDet==4){
          layer = tTopo->tidRing(detId);
          nStrips_tid[layer-1] += sgd->specificTopology().nstrips();
        } 
        else if(SubDet==5){
          layer = tTopo->tobLayer(detId);
          nStrips_tob[layer-1] += sgd->specificTopology().nstrips();
        } 
        else if(SubDet==6){
          layer = tTopo->tecRing(detId);
          nStrips_tec[layer-1] += sgd->specificTopology().nstrips();
        }


        const SiPileUpSignals::SignalMapType* theSignal(theSiPileUpSignals->getSignal(sgd->geographicalId().rawId()));  

        if(theSignal) {
          for(const auto& amp : *theSignal) {
            unsigned int minAmplitude = 10000;
            if ( amp.second > minAmplitude ) {
              if(SubDet==3) {
                ++nStripsNonZero_tib[layer-1];
                h_charge_tib.Fill( amp.second );
              }
              else if(SubDet==4){
                ++nStripsNonZero_tid[layer-1];
                h_charge_tid.Fill( amp.second );
              } 
              else if(SubDet==5){
                ++nStripsNonZero_tob[layer-1];
                h_charge_tob.Fill( amp.second );
              } 
              else if(SubDet==6){
                ++nStripsNonZero_tec[layer-1];
                h_charge_tec.Fill( amp.second );
              } 
            }
          }
        }
      }
    }


    for ( unsigned int i = 0; i < nStrips_tib.size(); ++i ) {
      occupancy_tib[i] = float(nStripsNonZero_tib[i]) / float(nStrips_tib[i]);
      generateAPVBaseline( occupancy_tib[i], h_charge_tib, apvBaselineDistributions_tib_[i] );
    }

    for ( unsigned int i = 0; i < nStrips_tob.size(); ++i ) {
      occupancy_tob[i] = float(nStripsNonZero_tob[i]) / float(nStrips_tob[i]);
      generateAPVBaseline( occupancy_tob[i], h_charge_tob, apvBaselineDistributions_tob_[i] );
    }

    for ( unsigned int i = 0; i < nStrips_tid.size(); ++i ) {
      occupancy_tid[i] = float(nStripsNonZero_tid[i]) / float(nStrips_tid[i]);
      generateAPVBaseline( occupancy_tid[i], h_charge_tid, apvBaselineDistributions_tid_[i] );
    }


    for ( unsigned int i = 0; i < nStrips_tec.size(); ++i ) {
      occupancy_tec[i] = float(nStripsNonZero_tec[i]) / float(nStrips_tec[i]);
      generateAPVBaseline( occupancy_tec[i], h_charge_tec, apvBaselineDistributions_tec_[i] );
    }
}

//============================================================================                

void SiStripDigitizerAlgorithm::generateAPVBaseline(
              float occupancy,
              TH1F chargeDistribution,
              TH1F& baselineDistribution
              ) {


    if ( occupancy == 0. ) {
      baselineDistribution.Fill(0);
      return;
    }

    unsigned int maxBXWithoutInteraction = log(0.01)/log(1-occupancy);  // Consider BX where probablitiy of interaction is > 1% of occupancy
    TF1 *f_probBX = new TF1("f_probBX","[0]*(1-[0])^x",0,maxBXWithoutInteraction);
    f_probBX->SetParameter(0,occupancy);

    for ( unsigned int n = 0; n < apv_nBaselineToGenerate; ++n ) {

      float baselineQ = 0;
      float timeSinceInteractionInMicroS = 0;
      unsigned int totalBX = 0;
      unsigned int totalInteractions = 0;
      unsigned int nInteractionsWithSmalChange = 0;

      for ( unsigned int i_interaction = 0; i_interaction < apv_nPreviousInteractionsToSimulate; ++i_interaction ) {
        ++totalInteractions;
        unsigned int BX = f_probBX->GetRandom();
        float charge = chargeDistribution.GetRandom();

        totalBX += BX;
        timeSinceInteractionInMicroS += float(totalBX) * 25 / 1000;

        float extraChargeToBaseline = ( charge * exp( -1.0 * timeSinceInteractionInMicroS / apv_decayConstantInMicroS ) );
        if ( extraChargeToBaseline / baselineQ < apv_smallChangeThreshold ) ++nInteractionsWithSmalChange; 
        
        if ( nInteractionsWithSmalChange > apv_smallChangeN ) break;

        baselineQ += extraChargeToBaseline;
      }

      float baselineV = baselineQ * apv_mVPerQ * apv_fCPerElectron;
      if ( baselineV > apv_maxResponse ) baselineV = apv_maxResponse;
      baselineDistribution.Fill( baselineV );
    }
}

//============================================================================                


void
SiStripDigitizerAlgorithm::digitize(
			   edm::DetSet<SiStripDigi>& outdigi,
			   edm::DetSet<SiStripRawDigi>& outrawdigi,
         edm::DetSet<SiStripRawDigi>& outStripAmplitudes,
         edm::DetSet<SiStripRawDigi>& outStripAmplitudesPostAPV,
         edm::DetSet<SiStripRawDigi>& outStripAPVBaselines,
			   edm::DetSet<StripDigiSimLink>& outLink,
			   const StripGeomDetUnit *det,
			   edm::ESHandle<SiStripGain> & gainHandle,
			   edm::ESHandle<SiStripThreshold> & thresholdHandle,
			   edm::ESHandle<SiStripNoises> & noiseHandle,
			   edm::ESHandle<SiStripPedestals> & pedestalHandle,
			   std::vector<std::pair<int,std::bitset<6>>> & theAffectedAPVvector,
                           CLHEP::HepRandomEngine* engine,
         const TrackerTopology *tTopo
          ) {
  unsigned int detID = det->geographicalId().rawId();
  int numStrips = (det->specificTopology()).nstrips();  

  // Moved from later (if CommonModeNoise)
  DetId  detId(detID);
  uint32_t SubDet = detId.subdetId();

  TH1F* apvBaselineDistribution = 0;
  int layer = -1;
  if(SubDet==3) {
    layer = tTopo->tibLayer(detId);
    apvBaselineDistribution = &apvBaselineDistributions_tib_[layer-1];
  }
  else if(SubDet==4){
    layer = tTopo->tidRing(detId);
    apvBaselineDistribution = &apvBaselineDistributions_tid_[layer-1];
  } 
  else if(SubDet==5){
    layer = tTopo->tobLayer(detId);
    apvBaselineDistribution = &apvBaselineDistributions_tob_[layer-1];
  } 
  else if(SubDet==6){
    layer = tTopo->tecRing(detId);
    apvBaselineDistribution = &apvBaselineDistributions_tec_[layer-1];
  } 


  const SiPileUpSignals::SignalMapType* theSignal(theSiPileUpSignals->getSignal(detID));  

  std::vector<float> detAmpl(numStrips, 0.);
  if(theSignal) {
    for(const auto& amp : *theSignal) {
      detAmpl[amp.first] = amp.second;
    }
  }

  //removing signal from the dead (and HIP effected) strips
  std::vector<bool>& badChannels = allBadChannels[detID];
  for(int strip =0; strip < numStrips; ++strip) {
    if(badChannels[strip]) {detAmpl[strip] = 0.;}
  }

  // SCD, before APV sim
  for(int strip =0; strip < numStrips; ++strip) {
    outStripAmplitudes.push_back(SiStripRawDigi(detAmpl[strip]/theElectronPerADC));;
  }

  // APV simulation here
  for(int strip =0; strip < numStrips; ++strip) {
    if (detAmpl[strip] > 0 ) {
      // Charge from electrons to fC
      double stripCharge = detAmpl[strip]*apv_fCPerElectron;

      double baselineV = apvBaselineDistribution->GetRandom();
      outStripAPVBaselines.push_back(SiStripRawDigi(baselineV));

      // Fitted parameters from G Hall/M Raymond
      double maxResponse = apv_maxResponse;
      double rate = apv_rate;
      double baselineQ = apv_maxResponse;

      // Convert V0 into baseline charge
      if ( baselineV < baselineQ ) {
        baselineQ = - 1.0 * rate * log( 2 * maxResponse / ( baselineV + maxResponse ) - 1);
      }

      // Add charge deposited in this BX
      double newStripCharge = baselineQ + stripCharge;

      // Apply APV response
      double signalV = 2 * maxResponse / ( 1 + exp( -1.0 * newStripCharge / rate) ) - maxResponse;
      double gain = signalV - baselineV;

      // Convert gain (mV) to charge (assuming linear region of APV) and then to electrons
      double outputCharge = gain/apv_mVPerQ;
      double outputChargeInADC = outputCharge / apv_fCPerElectron;

      detAmpl[strip] = outputChargeInADC;
    }
  }

  // SCD, after APV sim
  for(int strip =0; strip < numStrips; ++strip) outStripAmplitudesPostAPV.push_back(SiStripRawDigi(detAmpl[strip]/theElectronPerADC));;



  if(APVSaturationFromHIP){
    //Implementation of the proper charge scaling function. Need consider resaturation effect:
    //The probability map gives  the probability that at least one HIP happened during the last N bunch crossings (cfr APV recovery time).
    //The impact on the charge depends on the clostest HIP occurance (in terms of bunch crossing).
    //The function discribing the APV recovery is therefore the weighted average function which takes into account all possibilities of HIP occurances across the last bx's.

    // do this step here because we now have access to luminosity information
    if(FirstDigitize_) {

      for(std::map<int,float>::iterator iter = mapOfAPVprobabilities.begin(); iter != mapOfAPVprobabilities.end(); ++iter){
        std::bitset<6> bs;
        for(int Napv=0;Napv<6;Napv++){
          float cursor=CLHEP::RandFlat::shoot(engine);
          bs[Napv]=cursor < iter->second*APVSaturationProb_ ? 1:0;  //APVSaturationProb has been scaled by PU luminosity
        }
        SiStripTrackerAffectedAPVMap[iter->first]=bs;
      }

      NumberOfBxBetweenHIPandEvent=1e3;
      bool HasAtleastOneAffectedAPV=false;
      while(!HasAtleastOneAffectedAPV){
        for(int bx=floor(300.0/25.0);bx>0;bx--){ //Reminder: make these numbers not hard coded!!
          float temp=CLHEP::RandFlat::shoot(engine)<0.5?1:0;
          if(temp==1 && bx<NumberOfBxBetweenHIPandEvent){
            NumberOfBxBetweenHIPandEvent=bx;
            HasAtleastOneAffectedAPV=true;
          }
        }
      }

      FirstDigitize_ = false;
    }

    std::bitset<6> & bs=SiStripTrackerAffectedAPVMap[detID];

    if(bs.any()){
      // store this information so it can be saved to the event later
      theAffectedAPVvector.push_back(std::make_pair(detID,bs));

      if(!PreMixing_) {

    // Here below is the scaling function which describes the evolution of the baseline (i.e. how the charge is suppressed).
    // This must be replaced as soon as we have a proper modeling of the baseline evolution from VR runs
    float Shift=1-NumberOfBxBetweenHIPandEvent/floor(300.0/25.0); //Reminder: make these numbers not hardcoded!! 
    float randomX=CLHEP::RandFlat::shoot(engine);
    float scalingValue=(randomX-Shift)*10.0/7.0-3.0/7.0;

    for(int strip =0; strip < numStrips; ++strip) {
      if(!badChannels[strip] &&  bs[strip/128]==1){
        detAmpl[strip] *=scalingValue>0?scalingValue:0.0;
      }
    }
        }
      }
  }


  SiStripNoises::Range detNoiseRange = noiseHandle->getRange(detID);
  SiStripApvGain::Range detGainRange = gainHandle->getRange(detID);
  SiStripPedestals::Range detPedestalRange = pedestalHandle->getRange(detID);

// -----------------------------------------------------------

  auto& firstChannelWithSignal = firstChannelsWithSignal[detID];
  auto& lastChannelWithSignal = lastChannelsWithSignal[detID];
  auto iAssociationInfoByChannel=associationInfoForDetId_.find(detID); // Use an iterator so that I can easily remove it once finished

  if(zeroSuppression){

    //Adding the strip noise
    //------------------------------------------------------  
    if(noise){ 
                         
      if(SingleStripNoise){
//      std::cout<<"In SSN, detId="<<detID<<std::endl;
        std::vector<float> noiseRMSv; 
        noiseRMSv.clear(); 
        noiseRMSv.insert(noiseRMSv.begin(),numStrips,0.); 
        for(int strip=0; strip< numStrips; ++strip){ 
          if(!badChannels[strip]){
            float gainValue = gainHandle->getStripGain(strip, detGainRange); 
            noiseRMSv[strip] = (noiseHandle->getNoise(strip,detNoiseRange))* theElectronPerADC/gainValue;
            //std::cout<<"<SiStripDigitizerAlgorithm::digitize>: gainValue: "<<gainValue<<"\tnoiseRMSv["<<strip<<"]: "<<noiseRMSv[strip]<<std::endl;
          }
        }
         theSiNoiseAdder->addNoiseVR(detAmpl, noiseRMSv, engine);
      } else {
          int RefStrip = int(numStrips/2.);
          while(RefStrip<numStrips&&badChannels[RefStrip]){ //if the refstrip is bad, I move up to when I don't find it 
            RefStrip++;
          } 
          if(RefStrip<numStrips){
            float RefgainValue = gainHandle->getStripGain(RefStrip, detGainRange);
            float RefnoiseRMS = noiseHandle->getNoise(RefStrip,detNoiseRange) *theElectronPerADC/RefgainValue; 
          
            theSiNoiseAdder->addNoise(detAmpl,firstChannelWithSignal,lastChannelWithSignal,numStrips,RefnoiseRMS, engine);
            //std::cout<<"<SiStripDigitizerAlgorithm::digitize>: RefgainValue: "<<RefgainValue<<"\tRefnoiseRMS: "<<RefnoiseRMS<<std::endl;
          }
      }
    }//if noise

    DigitalVecType digis;
    theSiZeroSuppress->suppress(theSiDigitalConverter->convert(detAmpl, gainHandle, detID), digis, detID,noiseHandle,thresholdHandle);
    // Now do the association to truth. Note that if truth association was turned off in the configuration this map
    // will be empty and the iterator will always equal associationInfoForDetId_.end().
    if( iAssociationInfoByChannel!=associationInfoForDetId_.end() ) { // make sure the readings for this DetID aren't completely from noise
      for( const auto& iDigi : digis ) {
        auto& associationInfoByChannel=iAssociationInfoByChannel->second;
        const std::vector<AssociationInfo>& associationInfo=associationInfoByChannel[iDigi.channel()];

        // Need to find the total from all sim hits, because this might not be the same as the total
        // digitised due to noise or whatever.
        float totalSimADC=0;
        for( const auto& iAssociationInfo : associationInfo ) totalSimADC+=iAssociationInfo.contributionToADC;
        // Now I know that I can loop again and create the links
        for( const auto& iAssociationInfo : associationInfo ) {
          // Note simHitGlobalIndex used to have +1 because TrackerHitAssociator (the only place I can find this value being used)
          // expected counting to start at 1, not 0.  Now changed.
          outLink.push_back( StripDigiSimLink( iDigi.channel(), iAssociationInfo.trackID, iAssociationInfo.simHitGlobalIndex, iAssociationInfo.tofBin, iAssociationInfo.eventID, iAssociationInfo.contributionToADC/totalSimADC ) );
        } // end of loop over associationInfo
      } // end of loop over the digis
    } // end of check that iAssociationInfoByChannel is a valid iterator
    outdigi.data = digis;
  }//if zeroSuppression

  if(!zeroSuppression){
    //if(noise){
      // the constant pedestal offset is needed because
      //   negative adc counts are not allowed in case
      //   Pedestal and CMN subtraction is performed.
      //   The pedestal value read from the conditions
      //   is pedValue and after the pedestal subtraction
      //   the baseline is zero. The Common Mode Noise
      //   is not subtracted from the negative adc counts
      //   channels. Adding pedOffset the baseline is set
      //   to pedOffset after pedestal subtraction and CMN
      //   is subtracted to all the channels since none of
      //   them has negative adc value. The pedOffset is
      //   treated as a constant component in the CMN
      //   estimation and subtracted as CMN.
      
         
		//calculating the charge deposited on each APV and subtracting the shift
		//------------------------------------------------------
		if(BaselineShift){
		   theSiNoiseAdder->addBaselineShift(detAmpl, badChannels);
		}
		
		//Adding the strip noise
		//------------------------------------------------------						 
		if(noise){
		    std::vector<float> noiseRMSv;
			noiseRMSv.clear();
		    noiseRMSv.insert(noiseRMSv.begin(),numStrips,0.);
			
		    if(SingleStripNoise){
			    for(int strip=0; strip< numStrips; ++strip){
			  		if(!badChannels[strip]) noiseRMSv[strip] = (noiseHandle->getNoise(strip,detNoiseRange))* theElectronPerADC;
			  	}
			
	    	} else {
			    int RefStrip = 0; //int(numStrips/2.);
		    	    while(RefStrip<numStrips&&badChannels[RefStrip]){ //if the refstrip is bad, I move up to when I don't find it
					RefStrip++;
				}
				if(RefStrip<numStrips){
					float noiseRMS = noiseHandle->getNoise(RefStrip,detNoiseRange) *theElectronPerADC;
					for(int strip=0; strip< numStrips; ++strip){
			       		if(!badChannels[strip]) noiseRMSv[strip] = noiseRMS;
			  		}
				}
			}
			
                    theSiNoiseAdder->addNoiseVR(detAmpl, noiseRMSv, engine);
		}			
		
		//adding the CMN
		//------------------------------------------------------
        if(CommonModeNoise){
		  float cmnRMS = 0.;
		  if(SubDet==3){
		    cmnRMS = cmnRMStib;
		  }else if(SubDet==4){
		    cmnRMS = cmnRMStid;
		  }else if(SubDet==5){
		    cmnRMS = cmnRMStob;
		  }else if(SubDet==6){
		    cmnRMS = cmnRMStec;
		  }
		  cmnRMS *= theElectronPerADC;
                  theSiNoiseAdder->addCMNoise(detAmpl, cmnRMS, badChannels, engine);
		}
		
        		
		//Adding the pedestals
		//------------------------------------------------------
		
		std::vector<float> vPeds;
		vPeds.clear();
		vPeds.insert(vPeds.begin(),numStrips,0.);
		
		if(RealPedestals){
		    for(int strip=0; strip< numStrips; ++strip){
			   if(!badChannels[strip]) vPeds[strip] = (pedestalHandle->getPed(strip,detPedestalRange)+pedOffset)* theElectronPerADC;
		    }
        } else {
		    for(int strip=0; strip< numStrips; ++strip){
			  if(!badChannels[strip]) vPeds[strip] = pedOffset* theElectronPerADC;
			}
		}
		
		theSiNoiseAdder->addPedestals(detAmpl, vPeds);	
		
		 
	//if(!RealPedestals&&!CommonModeNoise&&!noise&&!BaselineShift&&!APVSaturationFromHIP){
    //  edm::LogWarning("SiStripDigitizer")<<"You are running the digitizer without Noise generation and without applying Zero Suppression. ARE YOU SURE???";
    //}else{							 
    
    DigitalRawVecType rawdigis = theSiDigitalConverter->convertRaw(detAmpl, gainHandle, detID);

    // Now do the association to truth. Note that if truth association was turned off in the configuration this map
    // will be empty and the iterator will always equal associationInfoForDetId_.end().
    if( iAssociationInfoByChannel!=associationInfoForDetId_.end() ) { // make sure the readings for this DetID aren't completely from noise
      // N.B. For the raw digis the channel is inferred from the position in the vector.
      // I'VE NOT TESTED THIS YET!!!!!
      // ToDo Test this properly.
      for( size_t channel=0; channel<rawdigis.size(); ++channel ) {
        auto& associationInfoByChannel=iAssociationInfoByChannel->second;
        const auto iAssociationInfo=associationInfoByChannel.find(channel);
        if( iAssociationInfo==associationInfoByChannel.end() ) continue; // Skip if there is no sim information for this channel (i.e. it's all noise)
        const std::vector<AssociationInfo>& associationInfo=iAssociationInfo->second;

        // Need to find the total from all sim hits, because this might not be the same as the total
        // digitised due to noise or whatever.
        float totalSimADC=0;
        for( const auto& iAssociationInfo : associationInfo ) totalSimADC+=iAssociationInfo.contributionToADC;
        // Now I know that I can loop again and create the links
        for( const auto& iAssociationInfo : associationInfo ) {
          // Note simHitGlobalIndex used to have +1 because TrackerHitAssociator (the only place I can find this value being used)
          // expected counting to start at 1, not 0.  Now changed.
          outLink.push_back( StripDigiSimLink( channel, iAssociationInfo.trackID, iAssociationInfo.simHitGlobalIndex, iAssociationInfo.tofBin, iAssociationInfo.eventID, iAssociationInfo.contributionToADC/totalSimADC ) );
        } // end of loop over associationInfo
      } // end of loop over the digis
    } // end of check that iAssociationInfoByChannel is a valid iterator

    outrawdigi.data = rawdigis;
	
	//}
  }

  // Now that I've finished with this entry in the map of associations, I can remove it.
  // Note that there might not be an association if the ADC reading is from noise in which
  // case associationIsValid will be false.
  if( iAssociationInfoByChannel!=associationInfoForDetId_.end() ) associationInfoForDetId_.erase(iAssociationInfoByChannel);
}
