#include "L1Trigger/L1TNtuples/interface/L1AnalysisPhaseIPFJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

L1Analysis::L1AnalysisPhaseIPFJet::L1AnalysisPhaseIPFJet()
{
}

L1Analysis::L1AnalysisPhaseIPFJet::~L1AnalysisPhaseIPFJet()
{

}



void L1Analysis::L1AnalysisPhaseIPFJet::SetPhaseIPFJet(const edm::Handle< vector<reco::CaloJet> > phaseIPFJets, unsigned maxL1Extra)
{

  for (unsigned int i=0; i<phaseIPFJets->size() && l1extra_.nPhaseIPFJets<maxL1Extra; i++){
    if (phaseIPFJets->at(i).pt()>0){
      if(phaseIPFJets->at(i).pt()>30 && abs(phaseIPFJets->at(i).eta())<2.4)
	   // l1extra_.phaseIPFJetHt += phaseIPFJets->at(i).pt();
      l1extra_.phaseIPFJetEt .push_back(phaseIPFJets->at(i).pt());
      l1extra_.phaseIPFJetEta.push_back(phaseIPFJets->at(i).eta());
      l1extra_.phaseIPFJetPhi.push_back(phaseIPFJets->at(i).phi());
      l1extra_.nPhaseIPFJets++;
    }
  }
}

void L1Analysis::L1AnalysisPhaseIPFJet::SetPhaseIPFSums  (const edm::Handle< BXVector<l1t::EtSum> >  phaseIPFSums)
{
  BXVector<l1t::EtSum>::const_iterator sumItr = phaseIPFSums->begin(0);
  const BXVector<l1t::EtSum>::const_iterator sumItrEnd = phaseIPFSums->end(0);
  for ( ; sumItr != sumItrEnd; ++sumItr ) {
    l1t::EtSum::EtSumType sumType{ sumItr->getType() };

    // std::cout << "Sum : " << sumItr->getType() << " " << sumItr->hwPt() << " " << sumItr->pt() << std::endl;

    if ( sumType == l1t::EtSum::kTotalHt ) {
      l1extra_.phaseIPFJetHt = sumItr->pt();
    }
    else if ( sumType == l1t::EtSum::kMissingEt ) {
      l1extra_.phaseIPFJetMET = sumItr->pt();
    }
  } 
}

void L1Analysis::L1AnalysisPhaseIPFJet::SetPFJet(const edm::Handle< vector<l1t::PFJet> > ak4PFJets, unsigned maxL1Extra)
{

  for (unsigned int i=0; i<ak4PFJets->size() && l1extra_.nPhaseIPFJets<maxL1Extra; i++){
    if (ak4PFJets->at(i).pt()>0){
      if(ak4PFJets->at(i).pt()>30 && abs(ak4PFJets->at(i).eta())<2.4)
	    l1extra_.ak4PFJetHt += ak4PFJets->at(i).pt();
      l1extra_.ak4PFJetEt .push_back(ak4PFJets->at(i).pt());
      l1extra_.ak4PFJetEta.push_back(ak4PFJets->at(i).eta());
      l1extra_.ak4PFJetPhi.push_back(ak4PFJets->at(i).phi());
      l1extra_.nAK4PFJets++;
    }
  }
}

void L1Analysis::L1AnalysisPhaseIPFJet::SetPFMET(const edm::Handle<reco::PFMETCollection> pfMET)
{
  l1extra_.pfMet = pfMET->at(0).pt();
}


void L1Analysis::L1AnalysisPhaseIPFJet::SetGenJet(const edm::Handle<reco::GenJetCollection> genJets, unsigned maxL1Extra)
{

  reco::GenJetCollection::const_iterator genJetItr = genJets->begin();
  reco::GenJetCollection::const_iterator genJetEnd = genJets->end();
  for( ; genJetItr != genJetEnd ; ++genJetItr) {
    if(genJetItr->pt()>30 && abs(genJetItr->eta())<2.4)
      l1extra_.genJetHt += genJetItr->pt();
    l1extra_.genJetPt.push_back( genJetItr->pt() );
    l1extra_.genJetEta.push_back( genJetItr->eta() );
    l1extra_.genJetPhi.push_back( genJetItr->phi() );
    l1extra_.genJetM.push_back( genJetItr->mass() );
    l1extra_.nGenJet++;
  }
}

void L1Analysis::L1AnalysisPhaseIPFJet::SetGenMET(const edm::Handle<reco::GenMETCollection> genMET)
{

  l1extra_.genMet = genMET->at(0).pt();
}


