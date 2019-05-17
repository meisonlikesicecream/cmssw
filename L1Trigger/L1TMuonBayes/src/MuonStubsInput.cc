/*
 * MuonStubsInput.cc
 *
 *  Created on: Jan 31, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#include "L1Trigger/L1TMuonBayes/interface/MuonStubsInput.h"

MuonStubsInput::MuonStubsInput(const ProcConfigurationBase* config): config(config), muonStubsInLayers(config->nLayers()) {

}


std::ostream & operator<< (std::ostream &out, const MuonStubsInput& stubsInput) {
  out <<"MuonStubsInput: "<<std::endl;
  for(auto& layerStubs : stubsInput.getMuonStubs()) {
    for(auto& stub : layerStubs) {
      out <<(*stub)<<std::endl;
    }
  }
  return out;
}



//gives stub phiHw or phiBHw - depending which layer is requested
//assumes that the banidg layer
const int MuonStubsInput::getPhiHw(unsigned int iLayer, unsigned int iInput) const {
  if(config->isBendingLayer(iLayer) ) {
    if(iInput >= muonStubsInLayers[iLayer -1].size())
      return MuonStub::EMTPY_PHI;
    return  muonStubsInLayers[iLayer -1][iInput]->phiBHw;
  }
  else {
    if(iInput >= muonStubsInLayers[iLayer].size())
      return MuonStub::EMTPY_PHI;
    return  muonStubsInLayers[iLayer][iInput]->phiBHw;
  }
}
