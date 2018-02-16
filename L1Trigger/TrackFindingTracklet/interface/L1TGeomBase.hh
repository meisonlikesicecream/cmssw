#ifndef L1TGEOMBASE_H
#define L1TGEOMBASE_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>
using namespace std;

#include "L1TStub.hh"
#include "L1TTracklet.hh"

#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/KF4ParamsComb.h"
#include "TMTrackTrigger/TMTrackFinder/interface/SimpleLR.h"


class L1TGeomBase{

  friend class L1TDisk;
  friend class L1TBarrel;

private:
  L1TGeomBase(){
  }


public:
  typedef std::map< L1TStub, edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ >  >, L1TStubCompare > stubMapType;

  void fitTracks() {
    for(int iSector=0;iSector<NSector_;iSector++){
      for(unsigned int i=0;i<tracklets_[iSector].size();i++){ 
	if (tracklets_[iSector][i].nStubs()>3){
	  L1TTrack aTrack(tracklets_[iSector][i]);
    tracks_[iSector].addTrack(aTrack);
	}
      }
    }
  }

  // Fit tracks using a method from TMTT
  void fitTracks( TMTT::Settings* settings, stubMapType stubMap, const TrackerGeometry* theTrackerGeom, const TrackerTopology* tTopo ) {

    TMTT::SimpleLR fitter(settings);
    fitter.initRun();
    // To run with the KF fitter, replace the previous two lines with the following
    // TMTT::KF4ParamsComb fitter(settings, 4, "KF4ParamsComb" );

    for(int iSector=0;iSector<NSector_;iSector++){
      for(unsigned int i=0;i<tracklets_[iSector].size();i++){ 
        if (tracklets_[iSector][i].nStubs()>3){
            // Convert the L1TStub on the L1TTracklet to TMTT::Stub
            L1TTracklet tracklet = tracklets_[iSector][i];
            vector<L1TStub> trackletStubs = tracklet.getStubs();
            vector<TMTT::Stub> tmttStubs;

            // Store map between L1TStub and TMTT::Stub
            // Useful for conversion back to tracklet class after fitting
            map<const TMTT::Stub*,const L1TStub*> trackletToTMTTStubMap;

            // Loop L1TStubs, and convert to TMTT:Stub
            stubMapType::const_iterator it;
            for (vector<L1TStub>::const_iterator itstubs = trackletStubs.begin(); itstubs != trackletStubs.end(); itstubs++) {
                it=stubMap.find(*itstubs);
                if (it!=stubMap.end()) {
                  TMTT::Stub tmttStub( it->second, 0, settings, theTrackerGeom, tTopo );
                  tmttStubs.push_back( tmttStub );
                  trackletToTMTTStubMap[&tmttStubs.back()] = &(*itstubs);
                }
            }

            // Convert vector of stubs to vector of pointers to stubs
            // Can probably be done in previous step
            vector<const TMTT::Stub*> tmttStubs_pointers;
            for (unsigned int j=0; j<tmttStubs.size(); ++j ){
              tmttStubs_pointers.push_back( &tmttStubs[j] );
            }

            // Gather info needed to create a TMTT::L1track3D
            // Store dummy values for HT cell location
            pair<unsigned int, unsigned int> cellLocationRphi(0,0);
            pair<unsigned int, unsigned int> cellLocationRz(0,0);
            // Helix parameters of L1track3D from tracklet
            pair<float, float> helixRphi(1/tracklet.pt(settings->getBfield()),tracklet.phi0());
            pair<float, float> helixRz(tracklet.z0(),tracklet.t());
            // Phi sector.  TMTT uses 8 (or 9), tracklet 24.
            // Assuming that using 24 phi sectors in a fitter from TMTT is ok
            unsigned int iPhiSec = iSector;
            // Eta regions
            // TMTT  uses 18 eta regions
            // The KF fitter needs to know which TMTT eta region the track is in
            // As it uses this information to know which layers the stubs are in (4 possibilities)
            // SimpleLR doesn't use this information, other than for debugging and for storing in output L1fittedTrack

            unsigned int iEtaReg = 8;
            double trackletEta = fabs(-log(tan(0.5*(0.25*8*atan(1.0)-atan(tracklet.t())))));
            if ( trackletEta > 2.16 ) iEtaReg = 0;
            else if ( trackletEta > 1.95 ) iEtaReg = 1;
            else if ( trackletEta > 1.7 ) iEtaReg = 2;
            else if ( trackletEta > 1.43 ) iEtaReg = 3;
            else if ( trackletEta > 1.16 ) iEtaReg = 4;
            else if ( trackletEta > 0.89 ) iEtaReg = 5;
            else if ( trackletEta > 0.61 ) iEtaReg = 6;
            else if ( trackletEta > 0.31 ) iEtaReg = 7;
            else iEtaReg = 8;

            // Create L1track3D object
            TMTT::L1track3D tmttTrack( settings, tmttStubs_pointers, cellLocationRphi, helixRphi, cellLocationRz, helixRz, iPhiSec, iEtaReg);

            // Fit track
            TMTT::L1fittedTrack tmttFittedTrack = fitter.fit(tmttTrack, iPhiSec, iEtaReg );

            // Convert L1fittedTrack to L1TTrack
            // Start by convert TMTT::Stubs back to tracklet L1Stub
            vector<L1TStub> l1TStubsOnFittedTrack;
            for (unsigned int j=0; j<tmttFittedTrack.getStubs().size(); ++j ){
              const TMTT::Stub* stub = tmttFittedTrack.getStubs()[j];

              if ( trackletToTMTTStubMap.find(stub) == trackletToTMTTStubMap.end() ) {
                // Some stubs appear to have a match in the map (same r, z coordinate), 
                // but have been modified (pointer to a different stub object)
                // Only noticed this happening with KF fitter, not explicitly checked with SimpleLR
                for (map<const TMTT::Stub*,const L1TStub*>::iterator i=trackletToTMTTStubMap.begin(); i!=trackletToTMTTStubMap.end(); ++i) {
                  double rDiff = fabs( i->first->r() - stub->r() );
                  double zDiff = fabs( i->first->z() - stub->z() );
                  if ( rDiff < 0.01 && zDiff < 0.01 ) {
                    l1TStubsOnFittedTrack.push_back( *i->second );
                  }
                }
              }
              else {
                l1TStubsOnFittedTrack.push_back( *trackletToTMTTStubMap[stub] );
              }

            }

            // Create L1TTrack object from L1fittedTrack and stubs
            // Set chi2 of fit to zero for now, to avoid any cuts later on
            L1TTrack aTrack(tracklet,
              tmttFittedTrack.qOverPt() * (0.00299792*settings->getBfield()),
              tmttFittedTrack.phi0(),
              tmttFittedTrack.z0(),
              tmttFittedTrack.tanLambda(),
              l1TStubsOnFittedTrack,
              0
              );

            tracks_[iSector].addTrack(aTrack);
        }
      }
    }
  }

  unsigned int nTracks(int iSector) const {return tracks_[iSector].size();}
  L1TTrack& track(int iSector, unsigned int i) {return tracks_[iSector].get(i);}

  unsigned int nTracklets(int iSector) const {return tracklets_[iSector].size();}
  L1TTracklet& tracklet(int iSector, unsigned int i) {return tracklets_[iSector][i];}

  L1TTracks allTracks(){
    L1TTracks tracks=tracks_[0];
    for (int i=1;i<NSector_;i++){
      tracks.addTracks(tracks_[i]);
    }
    return tracks;
  }

private:

  int NSector_;

  vector<vector<L1TStub> > stubs_;
  vector<vector<L1TTracklet> > tracklets_;
  vector<L1TTracks > tracks_;

};

#endif

