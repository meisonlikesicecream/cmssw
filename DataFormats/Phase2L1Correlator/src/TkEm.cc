// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     TkEm
//

#include "DataFormats/Phase2L1Correlator/interface/TkEm.h"

using namespace l1t;

TkEm::TkEm() {}

TkEm::TkEm(const LorentzVector& p4, const edm::Ref<EGammaBxCollection>& egRef, float tkisol)
    : L1Candidate(p4), egRef_(egRef), TrkIsol_(tkisol), TrkIsolPV_(-999) {}

TkEm::TkEm(const LorentzVector& p4,
                               const edm::Ref<EGammaBxCollection>& egRef,
                               float tkisol,
                               float tkisolPV)
    : L1Candidate(p4), egRef_(egRef), TrkIsol_(tkisol), TrkIsolPV_(tkisolPV) {}
