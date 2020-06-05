#include <iostream>
#include <string>
#include <vector>

#include "TTree.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "Geometry/HGCalCommonData/interface/HGCalGeometryMode.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"

#include "Geometry/HGCalCommonData/interface/HGCalWaferIndex.h"


#include <cstdlib>

namespace {
  template <typename T>
  struct array_deleter {
    void operator()(T* arr) { delete[] arr; }
  };
}  // namespace

class HGCalGeomTester : public edm::stream::EDAnalyzer<> {
public:
  explicit HGCalGeomTester(const edm::ParameterSet&);
  ~HGCalGeomTester();

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

private:
  void fillGeometry( const edm::EventSetup& es );
  void fillGeometry( const HGCalGeometry* geometry );

  void fillTriggerGeometry( const edm::EventSetup& es );
  void fillTriggerGeometry( const HGCalGeometry* geometry, const HGCalTopology& topology );

  edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> geomTokenEE_;
  edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> geomTokenHE_;

  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;

  edm::Service<TFileService> fs_;

  TTree* treeTriggerCells_;
  TTree* treeCells_;

  // tree variables
  int triggerCellSide_ = 0;
  int triggerCellSubdet_ = 0;
  int triggerCellLayer_ = 0;
  int triggerCellWaferU_ = 0;
  int triggerCellWaferV_ = 0;
  int triggerCellWaferType_ = 0;
  bool triggerCellValidDetId_ = 0;
  int triggerCellWaferShape_ = 0;
  int triggerCellWaferShapeRotation_ = 0;
  int triggerCellU_ = 0;
  int triggerCellV_ = 0;

  int cellSide_ = 0;
  int cellSubdet_ = 0;
  int cellLayer_ = 0;
  int cellWaferU_ = 0;
  int cellWaferV_ = 0;
  int cellWaferType_ = 0;
  bool cellValidDetId_ = false;
  int cellWaferShape_ = 0;
  int cellWaferShapeRotation_ = 0;
  int cellU_ = 0;
  int cellV_ = 0;

private:
  typedef std::unordered_map<uint32_t, std::unordered_set<uint32_t>> trigger_map_set;
};

/*****************************************************************/
HGCalGeomTester::HGCalGeomTester(const edm::ParameterSet& conf)
/*****************************************************************/
    : geomTokenEE_{esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", "HGCalEESensitive"})},
      geomTokenHE_{esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", "HGCalHESiliconSensitive"})}
{
  // // initialize output trees
  // //
  treeTriggerCells_ = fs_->make<TTree>("TreeTriggerCells", "Tree of all HGC trigger cells");
  // treeTriggerCells_->Branch("id", &triggerCellId_, "id/I");
  treeTriggerCells_->Branch("zside", &triggerCellSide_, "zside/I");
  treeTriggerCells_->Branch("subdet", &triggerCellSubdet_, "subdet/I");
  treeTriggerCells_->Branch("layer", &triggerCellLayer_, "layer/I");
  treeTriggerCells_->Branch("waferu", &triggerCellWaferU_, "waferu/I");
  treeTriggerCells_->Branch("waferv", &triggerCellWaferV_, "waferv/I");
  treeTriggerCells_->Branch("wafertype", &triggerCellWaferType_, "wafertype/I");
  treeTriggerCells_->Branch("validDetId", &triggerCellValidDetId_, "validDetId/I");
  treeTriggerCells_->Branch("wafershape", &triggerCellWaferShape_, "wafershape/I");
  treeTriggerCells_->Branch("wafershaperotation", &triggerCellWaferShapeRotation_, "wafershaperotation/I");
  treeTriggerCells_->Branch("triggercellu", &triggerCellU_, "triggercellu/I");
  treeTriggerCells_->Branch("triggercellv", &triggerCellV_, "triggercellv/I");
  // //
  treeCells_ = fs_->make<TTree>("TreeCells", "Tree of all HGC cells");
  treeCells_->Branch("zside", &cellSide_, "zside/I");
  treeCells_->Branch("subdet", &cellSubdet_, "subdet/I");
  treeCells_->Branch("layer", &cellLayer_, "layer/I");
  treeCells_->Branch("waferu", &cellWaferU_, "waferu/I");
  treeCells_->Branch("waferv", &cellWaferV_, "waferv/I");
  treeCells_->Branch("wafertype", &cellWaferType_, "wafertype/I");
  treeCells_->Branch("validDetId", &cellValidDetId_, "validDetId/I");
  treeCells_->Branch("wafershape", &cellWaferShape_, "wafershape/I");
  treeCells_->Branch("wafershaperotation", &cellWaferShapeRotation_, "wafershaperotation/I");
  treeCells_->Branch("cellu", &cellU_, "cellu/I");
  treeCells_->Branch("cellv", &cellV_, "cellv/I");
  // //
}

/*****************************************************************/
HGCalGeomTester::~HGCalGeomTester()
/*****************************************************************/
{}

/*****************************************************************/
void HGCalGeomTester::beginRun(const edm::Run& /*run*/, const edm::EventSetup& es)
/*****************************************************************/
{
}

/*****************************************************************/
void HGCalGeomTester::fillGeometry( const edm::EventSetup& es )
/*****************************************************************/
{
    // es.get<IdealGeometryRecord>().get("", geometry_);
  const auto& geomR_EE = es.getData(geomTokenEE_);
  const HGCalGeometry* geometry_EE = &geomR_EE;

  fillGeometry( geometry_EE );


  const auto& geomR_HE = es.getData(geomTokenHE_);
  const HGCalGeometry* geometry_HE = &geomR_HE;

  fillGeometry( geometry_HE );
}

/*****************************************************************/
void HGCalGeomTester::fillGeometry( const HGCalGeometry* geometry )
/*****************************************************************/
{
  const std::vector<DetId>& ids = geometry->getValidDetIds();

  for (const auto& id : ids ) {
    HGCSiliconDetId detid(id);

    cellSide_ = detid.zside();
    cellSubdet_ = detid.subdet();
    cellLayer_ = detid.layer();
    cellWaferU_ = detid.waferU();
    cellWaferV_ = detid.waferV();
    cellU_ = detid.cellU();
    cellV_ = detid.cellV();
    cellWaferType_ = detid.type();
    cellValidDetId_ = geometry->topology().valid(id);

    std::pair<int, int> typeAndRot{ geometry->topology().dddConstants().waferTypeRotation( cellLayer_, cellWaferU_, cellWaferV_ ) };
    cellWaferShape_ = typeAndRot.first;
    cellWaferShapeRotation_ = typeAndRot.second;

    // Print information on cells in three modules in the first CE-H layer
    // The three modules should be identical, whose location differ by a rotation by 120 (or 240) degrees
    if ( cellLayer_ == 1 && cellSubdet_ == 9 && cellSide_ > 0 ) {
      if ( 
            ( cellWaferU_ == 0 && cellWaferV_ == 2 ) || // This wafer has cellWaferShape_ = 4, and cellValidDetId_ = False
            ( cellWaferU_ == -3 && cellWaferV_ == -3 ) || // This wafer has cellWaferShape_ = 99, and cellValidDetId_ = True
            ( cellWaferU_ == 2 && cellWaferV_ == -1 ) // This wafer has cellWaferShape_ = 99, and cellValidDetId_ = True
       ) {
        std::cout << "Cell info, subdet : " << cellSubdet_ << ", layer : " << cellLayer_ << ", wafer u : " << cellWaferU_ << ", wafer v : " << cellWaferV_  << ", wafer shape : " << cellWaferShape_ << ", is valid :  " << cellValidDetId_ << std::endl;
      }
    }

    treeCells_->Fill();
  }
}

/*****************************************************************/
void HGCalGeomTester::fillTriggerGeometry( const edm::EventSetup& es )
/*****************************************************************/
{
    es.get<CaloGeometryRecord>().get("", triggerGeometry_);

    fillTriggerGeometry( triggerGeometry_->eeGeometry(), triggerGeometry_->eeTopology() );
    fillTriggerGeometry( triggerGeometry_->hsiGeometry(), triggerGeometry_->hsiTopology() );
}

/*****************************************************************/
void HGCalGeomTester::fillTriggerGeometry( const HGCalGeometry* geometry, const HGCalTopology& topology )
/*****************************************************************/
{
  // const HGCalParameters* hgcalParams = topology.dddConstants().getParameter();

  for (const auto& id : geometry->getValidDetIds()) {
    HGCalTriggerDetId id_trig( triggerGeometry_->getTriggerCellFromCell(id) );
    triggerCellSide_ = id_trig.zside();
    triggerCellSubdet_ = id_trig.subdet();
    triggerCellLayer_ = id_trig.layer();
    triggerCellWaferU_ = id_trig.waferU();
    triggerCellWaferV_ = id_trig.waferV();
    triggerCellWaferType_ = id_trig.type();
    triggerCellValidDetId_ = topology.valid(id);
    triggerCellU_ = id_trig.triggerCellU();
    triggerCellV_ = id_trig.triggerCellV();

    std::pair<int, int> typeAndRot{ geometry->topology().dddConstants().waferTypeRotation( triggerCellLayer_, triggerCellWaferU_, triggerCellWaferV_ ) };
    triggerCellWaferShape_ = typeAndRot.first;
    triggerCellWaferShapeRotation_ = typeAndRot.second;

    treeTriggerCells_->Fill();
  }
}

/*****************************************************************/
void HGCalGeomTester::analyze(const edm::Event& e, const edm::EventSetup& es)
/*****************************************************************/
{
  fillGeometry( es );

  // fillTriggerGeometry( es );
}

// define this as a plug-in
DEFINE_FWK_MODULE(HGCalGeomTester);
