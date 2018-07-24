#include "TMTrackTrigger/TMTrackFinder/interface/TrackerGeometryInfo.h"
#include <TGraph.h>
#include <TF1.h>

namespace TMTT {

TrackerGeometryInfo::TrackerGeometryInfo():
  barrelNTiltedModules_(12),
  barrelNLayersWithTiltedModules_(3),
  moduleZoR_(),
  moduleB_()
  {

  }

void TrackerGeometryInfo::getTiltedModuleInfo(const Settings* settings, const TrackerTopology*  theTrackerTopo, const TrackerGeometry* theTrackerGeom) {

  for (const GeomDetUnit* gd : theTrackerGeom->detUnits()) {
    float R0 = gd->position().perp();
    float Z0 = fabs(gd->position().z());
    DetId detId = gd->geographicalId();

    bool psModule = theTrackerGeom->getDetectorType( detId ) ==  TrackerGeometry::ModuleType::Ph2PSP;
    bool barrel = detId.subdetId()==StripSubdetector::TOB || detId.subdetId()==StripSubdetector::TIB;
    bool tiltedBarrel = barrel && (theTrackerTopo->tobSide(detId) != 3);
    bool isLower = theTrackerTopo->isLower(detId);

    // Calculate and store r, z, B at centre of titled modules
    // Only do this for the "lower" sensor
    if ( tiltedBarrel && isLower ) {

      int type   = 2*theTrackerTopo->tobSide(detId)-3; // -1 for tilted-, 1 for tilted+, 3 for flat
      int layer = theTrackerTopo->layer(detId);
      int ladder = theTrackerTopo->tobRod(detId);
      double corr   = (barrelNTiltedModules_+1)/2.;

      ladder = corr-(corr-ladder)*type; // Corrected ring number, bet 0 and barrelNTiltedModules_.at(layer), in ascending |z|

      double minZOfThisModule = fabs(gd->surface().zSpan().first);
      double maxZOfThisModule = fabs(gd->surface().zSpan().second);
      double minROfThisModule = gd->surface().rSpan().first;
      double maxROfThisModule = gd->surface().rSpan().second;

      // Calculate module tilt and B at centre of module
      double moduleTilt = atan( fabs(maxZOfThisModule - minZOfThisModule) / ( maxROfThisModule - minROfThisModule )  );
      double moduleTheta = atan(R0 / Z0);
      double BCorrection = fabs( cos( fabs(moduleTheta) - moduleTilt ) / sin( moduleTheta ) );

      // Only store if this value of B has not been stored already
      bool storeThisBCorrection = true;
      for (unsigned int iCorr = 0; iCorr < moduleB_.size(); ++iCorr ) {
        if ( fabs( BCorrection - moduleB_[iCorr] ) < 0.0001 ) {
          storeThisBCorrection = false;
          break;
        }
      }

      if ( storeThisBCorrection ) {
        moduleZoR_.push_back( Z0 / R0 );
        moduleB_.push_back( BCorrection );
      }
    }
  }

}

}
