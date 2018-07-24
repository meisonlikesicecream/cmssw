#include "TMTrackTrigger/TMTrackFinder/interface/KillDupFitTrks.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <map>

namespace TMTT {

//=== Make available cfg parameters & specify which algorithm is to be used for duplicate track removal.

void KillDupFitTrks::init(const Settings* settings, unsigned int dupTrkAlg)
{
  settings_ = settings;
  dupTrkAlg_ = dupTrkAlg;
  killDupTrks_.init(settings, dupTrkAlg); // Initialise duplicate removal algorithms that are common to all tracks.
}

//=== Eliminate duplicate tracks from the input collection, and so return a reduced list of tracks.

vector<L1fittedTrack> KillDupFitTrks::filter(const vector<L1fittedTrack>& vecTracks) const
{
  if (dupTrkAlg_ == 0) {

    // We are not running duplicate removal, so return original fitted track collection.    
    return vecTracks;

  } else {

    // Choose which algorithm to run, based on parameter dupTrkAlg_.
    switch (dupTrkAlg_) {
      // Run filters that only work on fitted tracks.
      case 50:  return filterAlg50( vecTracks )  ; break;
      case 51:  return filterAlg51( vecTracks ); break;
      // Run filters that work on any type of track (l1track2d, l1track3d, l1fittedtrack). 
      default:  return killDupTrks_.filter(vecTracks); 
    }

    // Should never get here ...	
  }
}

//=== Duplicate removal algorithm designed to run after the track helix fit, which eliminates duplicates  
//=== simply by requiring that the fitted (q/Pt, phi0) of the track correspond to the same HT cell in 
//=== which the track was originally found by the HT.
//=== N.B. This code runs on tracks in a single sector. It could be extended to run on tracks in entire
//=== tracker by adding the track's sector number to memory "htCellUsed" below.


vector<L1fittedTrack> KillDupFitTrks::filterAlg50(const vector<L1fittedTrack>& tracks) const
{
  // Hard-wired options to play with.
  const bool debug = false;
  const bool doRecoveryStep = true; // Do 2nd pass through rejected tracks to see if any should be rescued.
  const bool reduceDups = true; // Option attempting to reduce duplicate tracks during 2nd pass.
  const bool memorizeAllHTcells = false; // First pass stores in memory all cells that the HT found tracks in, not just those of tracks accepted by the first pass.
  const bool doSectorCheck = false; // Require fitted helix to lie within sector.
  const bool usePtAndZ0Cuts = false;
  const bool goOutsideArray = false; // Also store in memory stubs outside the HT array during 2nd pass.

  if (debug && tracks.size() > 0) cout<<"START "<<tracks.size()<<endl;

  vector<L1fittedTrack> tracksFiltered;

  // Make a first pass through the tracks, doing initial identification of duplicate tracks.
  set< pair<unsigned int, unsigned int> > htCellUsed;
  vector<const L1fittedTrack*> tracksRejected;

  for (const L1fittedTrack& trk : tracks) {
    // Only consider tracks whose fitted helix parameters are in the same sector as the HT originally used to find the track.
    if ( ( ! doSectorCheck) || trk.consistentSector() ) {
      if ( ( ! usePtAndZ0Cuts) || ( fabs(trk.z0()) < settings_->beamWindowZ() && trk.pt() > settings_->houghMinPt() - 0.2) ) {
	// Check if this track's fitted (q/pt, phi0) helix parameters correspond to the same HT cell as the HT originally found the track in.
	bool consistentCell = trk.consistentHTcell();
	if (consistentCell) {
	  // This track is probably not a duplicate, so keep & and store its HT cell location (which equals the HT cell corresponding to the fitted track).
	  tracksFiltered.push_back(trk);
	  // Memorize HT cell location corresponding to this track (identical for HT track & fitted track).
	  if ( ! memorizeAllHTcells) htCellUsed.insert( trk.getCellLocationHT() );
	  if (debug) {
	    const TP* tp = trk.getMatchedTP();
	    int tpIndex = (tp != nullptr) ? tp->index() : -999;
	    cout<<"FIRST PASS: m="<<trk.getCellLocationHT().first<<"/"<<trk.getCellLocationFit().first<<" c="<<trk.getCellLocationHT().second<<"/"<<trk.getCellLocationFit().second<<" tp="<<tpIndex<<" pure="<<trk.getPurity()<<" merged="<<trk.getL1track3D().mergedHTcell()<<endl;
	  }
	} else {
	  tracksRejected.push_back(&trk);
	}
	// Memorize HT cell location corresponding to this track, even if it was not accepted by first pass..
	if (memorizeAllHTcells) htCellUsed.insert( trk.getCellLocationHT() );
      }
    }
  }

  if (doRecoveryStep) {
    // Making a second pass through the rejected tracks, checking if any should be rescued.
    for (const L1fittedTrack* trk : tracksRejected) {
      // Get location in HT array corresponding to fitted track helix parameters.
      pair<unsigned int, unsigned int> htCell = trk->getCellLocationFit();
      // If this HT cell was not already memorized, rescue this track, since it is probably not a duplicate,
      // but just a track whose fitted helix parameters are a bit wierd for some reason.
      if (std::count(htCellUsed.begin(), htCellUsed.end(), htCell) == 0) {
	// If this track is from 2x2 merged HT cell, then check no track was found in any of 2x2 cells.
	bool mergedOK = true;
	if (trk->getL1track3D().mergedHTcell()) {
	  pair<unsigned int, unsigned int> htCell10( htCell.first + 1, htCell.second);
	  pair<unsigned int, unsigned int> htCell01( htCell.first    , htCell.second + 1);
	  pair<unsigned int, unsigned int> htCell11( htCell.first + 1, htCell.second + 1);
	  mergedOK = (std::count(htCellUsed.begin(), htCellUsed.end(), htCell10) == 0 &&
	              std::count(htCellUsed.begin(), htCellUsed.end(), htCell01) == 0 &&
	              std::count(htCellUsed.begin(), htCellUsed.end(), htCell11) == 0);
	}
	if (mergedOK) {
	  tracksFiltered.push_back(*trk); // Rescue track.
	  // Optionally store cell location to avoid rescuing other tracks at the same location, which may be duplicates of this track. 
	  bool outsideCheck =( goOutsideArray || trk->pt() > settings_->houghMinPt() );
	  if (reduceDups && outsideCheck) htCellUsed.insert( htCell );
	  if (debug) {
	    const TP* tp = trk->getMatchedTP();
	    int tpIndex = (tp != nullptr) ? tp->index() : -999;
	    cout<<"SECOND PASS: m="<<trk->getCellLocationHT().first<<"/"<<trk->getCellLocationFit().first<<" c="<<trk->getCellLocationHT().second<<"/"<<trk->getCellLocationFit().second<<" tp="<<tpIndex<<" pure="<<trk->getPurity()<<" merged="<<trk->getL1track3D().mergedHTcell()<<endl;
	  }
	}
      }
    }
  }

  // Debug printout to identify duplicate tracks that survived.
  if (debug) this->printDuplicateTracks(tracksFiltered);

  return tracksFiltered;
}
//=== Duplicate removal algorithm designed to run after the track helix fit, which eliminates duplicates  
//=== simply by requiring that no two tracks should have fitted (q/Pt, phi0) that correspond to the same HT
//=== cell. If they do, then only the first to arrive is kept.
//=== N.B. This code runs on tracks in a single sector. It could be extended to run on tracks in entire
//=== tracker by adding the track's sector number to memory "htCellUsed" below.

vector<L1fittedTrack> KillDupFitTrks::filterAlg51(const vector<L1fittedTrack>& tracks) const
{
  // Hard-wired options to play with.
  const bool debug = false;

  if (debug && tracks.size() > 0) cout<<"START "<<tracks.size()<<endl;

  vector<L1fittedTrack> tracksFiltered;
  set< pair<unsigned int, unsigned int> > htCellUsed;

  for (const L1fittedTrack& trk : tracks) {
      // Get location in HT array corresponding to fitted track helix parameters.
      pair<unsigned int, unsigned int> htCell = trk.getCellLocationFit();
      // If this HT cell was not already memorized, rescue this track, since it is probably not a duplicate,
      // but just a track whose fitted helix parameters are a bit wierd for some reason.
      if (std::count(htCellUsed.begin(), htCellUsed.end(), htCell) == 0) {
	tracksFiltered.push_back(trk); // Rescue track.
	// Store cell location to avoid rescuing other tracks at the same location, which may be duplicates of this track. 
	htCellUsed.insert( htCell );
	if (debug) {
	  const TP* tp = trk.getMatchedTP();
	  int tpIndex = (tp != nullptr) ? tp->index() : -999;
	  cout<<"ALG51: m="<<trk.getCellLocationHT().first<<"/"<<trk.getCellLocationFit().first<<" c="<<trk.getCellLocationHT().second<<"/"<<trk.getCellLocationFit().second<<" tp="<<tpIndex<<" pure="<<trk.getPurity()<<endl;
	}
      }
    }

  // Debug printout to identify duplicate tracks that survived.
  if (debug) this->printDuplicateTracks(tracksFiltered);

  return tracksFiltered;
}

// Debug printout of which tracks are duplicates.
void KillDupFitTrks::printDuplicateTracks(const vector<L1fittedTrack>& tracks) const {
  map<const TP*, vector<const L1fittedTrack*> > tpMap;
  for (const L1fittedTrack& trk : tracks) {
    const TP* tp = trk.getMatchedTP();
    if (tp != nullptr) {
      tpMap[tp].push_back(&trk);
    }
  }
  for (const auto& p : tpMap) {
    const TP* tp     = p.first;
    const vector<const L1fittedTrack*> vecTrk = p.second;
    if (vecTrk.size() > 1) {
      for (const L1fittedTrack* trk : vecTrk) {
	cout<<"  MESS UP : m="<<trk->getCellLocationHT().first<<"/"<<trk->getCellLocationFit().first<<" c="<<trk->getCellLocationHT().second<<"/"<<trk->getCellLocationFit().second<<" tp="<<tp->index()<<" pure="<<trk->getPurity()<<endl;
	cout<<"     stubs = ";
	for (const Stub* s : trk->getStubs()) cout<<s->index()<<" ";
	cout<<endl;
      }
    }
  }
  if (tracks.size() > 0) cout<<"FOUND "<<tracks.size()<<endl;
}

}
