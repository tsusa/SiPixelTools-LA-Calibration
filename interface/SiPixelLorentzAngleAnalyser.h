#ifndef CalibTracker_SiPixelLorentzAngle_SiPixelLorentzAngleAnalyser_h
#define CalibTracker_SiPixelLorentzAngle_SiPixelLorentzAngleAnalyser_h

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h" 
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
// #include "CalibTracker/SiPixelLorentzAngle/interface/TrackLocalAngle.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include <Geometry/CommonTopologies/interface/PixelGeomDetUnit.h>
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "CalibTracker/Records/interface/SiPixelTemplateDBObjectESProducerRcd.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelTemplateDBObject.h"
#include "CondFormats/SiPixelTransient/interface/SiPixelTemplate.h"
#include "CondFormats/SiPixelTransient/interface/SiPixelTemplateDefs.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TROOT.h>
#include "LAStructures.h"
#include "MCSHistograms.h"

/**
* 
* Class to determine the lorentz angle in the barrel pixel detector 
* returns a txt file with the fit for every of the 8 rings in the 3 layers
* 
*/

// ggiurgiu@fnal.gov : remove namespace 12/27/09
//namespace
//{
/*
  static const int maxpix = 1000;
  struct Pixinfo
  {
    int npix;
    float row[maxpix];
    float col[maxpix];
    float adc[maxpix];
    float x[maxpix];
    float y[maxpix];
  };
struct Hit
  {
    float x;
    float y;
    double alpha;
    double beta;
    double gamma;
  };
  struct Clust 
  {
    float x;
    float y;
    float charge;
    int size_x;
    int size_y;
    int maxPixelCol;
    int maxPixelRow;
    int minPixelCol;
    int minPixelRow;
  };
  struct Rechit 
  {
    float x;
    float y;
  };
*/
//}

class SiPixelLorentzAngleAnalyser : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> 
{
public:

  explicit SiPixelLorentzAngleAnalyser(const edm::ParameterSet& conf);

  virtual ~SiPixelLorentzAngleAnalyser();
  virtual void beginJob();
  virtual void endJob();
  void beginRun(const edm::Run& r, const edm::EventSetup& c) override;
  void endRun(const edm::Run& r, const edm::EventSetup& c) override {}  
  void analyze(const edm::Event& e, const edm::EventSetup& c) override;

 private:
  void ntupleAnalysisFpix();
  int fpixMinimumClusterSize();

  //int lower_bin_;  
  edm::Service<TFileService> fs;
  void fillPix(const SiPixelCluster & LocPix, const PixelTopology * topol, Pixinfo& pixinfo);
  void surface_derormation(const PixelTopology *topol, 
			   TrajectoryStateOnSurface &tsos, 
			   const SiPixelRecHit *recHitPix,
			   LocalPoint &lp_track, 
			   LocalPoint &lp_rechit);
  
  // tree branches barrel
  int run_;
  ULong64_t event_;
  int lumiblock_;
  int bx_;
  int orbit_;
  int module_;
  int ladder_;
  int layer_;
  int isflipped_;
  float pt_;
  float eta_;
  float phi_;
  double chi2_;
  double ndof_;
  Pixinfo pixinfo_;
  Hit simhit_, trackhit_;
  Clust clust_;
  Rechit rechit_;
  Rechit rechitCorr_;
  float trackhitCorrX_;
  float trackhitCorrY_;
  float qScale_;
  float rQmQt_;
  
  // tree branches forward
  int runF_;
  long int eventF_;  
  int sideF_;
  int diskF_;
  int bladeF_;
  int panelF_;
  int moduleF_;
  float ptF_;
  float etaF_;
  float phiF_;
  double chi2F_;
  double ndofF_;
  Pixinfo pixinfoF_;
  Hit simhitF_, trackhitF_;
  Clust clustF_;
  Rechit rechitF_;
  Rechit rechitCorrF_;
  float trackhitCorrXF_;
  float trackhitCorrYF_;
  float qScaleF_;
  float rQmQtF_;
  float magFieldF_[3];
  
  // parameters from config file
  int year_;
  std::string filename_;
  std::string filenameFit_;
  std::string source_;
  std::string subdetector_;
  std::string analysisType_;
  double fitRegionMin_;
  double fitRegionMax_;
  double fitRegionLinMin_;
  double fitRegionLinMax_;
  int nBinsX_;
  double xMin_, xMax_;
  int nBinsY_;
  double yMin_, yMax_;
  double ptmin_;
  double normChi2Max_;
  std::vector<int> clustSizeYMin_;
  int clustSizeYMax_;
  int clustSizeXMax_;
  double residualMax_;
  double clustChargeMax_;
  bool largePixCut_;
  int maxNumberOfPixelsCut_;  
  std::vector<std::string> root_fileNames_;

  const double width_ = 0.0285;
  const double ypitch_ = 0.0150;
  const double hypitch_ = ypitch_/2.;  
  
  // ==== Histograms
  // ----  Constants
  static const int nLayersMax_ = 4;
  static const int nModulesMax_ = 8;
  static const int nZSidesMax_ = 2;
  static const int nFlippedMax_ = 2;
  static const int nCentralEndMax_ = 2;
  // FPix 
  static const int nRingsMax_ = 2;
  static const int nPanelsMax_ = 2;
  
  // --- Control Distributions
  TH1F *h_ptBC_, *h_ptAC_;
  TH1F *h_ndofBC_, *h_ndofAC_;
  TH1F *h_chi2BC_, *h_chi2AC_;
  TH1F *h_normChi2BC_, *h_normChi2AC_;
  TH1F *h_nPixBC_, *h_nPixAC_;
  TH1F *h_cotAlphaBC_, *h_cotAlphaAC_;
  TH1F *h_cotBetaBC_, *h_cotBetaAC_;
  TH1F *h_cotBetaLayerBC_[nLayersMax_], *h_cotBetaLayerAC_[nLayersMax_];
  TH1F *h_largePixBC_, *h_largePixAC_;
  TH1F *h_clusterSizeYBC_, *h_clusterSizeYAC_;
  TH1F *h_clusterSizeXBC_, *h_clusterSizeXAC_;
  TH1F *h_clusterSizeYLayerBC_[nLayersMax_];
  TH1F *h_clusterSizeYLayerAC_[nLayersMax_];
  TH1F *h_clusterChargeBC_, *h_clusterChargeAC_;
  TH1F *h_residualBC_, *h_residualAC_;
  
  MCSHistograms mcshFPixAlpha_[nRingsMax_][nPanelsMax_][nZSidesMax_];
  MCSHistograms mcshFPixBeta_[nRingsMax_][nPanelsMax_][nZSidesMax_];
  TH1D  *magFieldHisto_[nRingsMax_][nPanelsMax_][nZSidesMax_][3];
  double meanMagField_[nRingsMax_][nPanelsMax_][nZSidesMax_][3];
  double countMagField_[nRingsMax_][nPanelsMax_][nZSidesMax_];
  
  TrackerHitAssociator::Config trackerHitAssociatorConfig_;
  
  MCSHistograms mcshTest_;
  int event_counter_, trackEventsCounter_,pixelTracksCounter_, hitCounter_, usedHitCounter_;
  
  int nlay;
  int nModules_[4];
  
  // CMSSW classes needed
  PropagatorWithMaterial  *thePropagator;
  PropagatorWithMaterial  *thePropagatorOp;
  KFUpdator *theUpdator;
  Chi2MeasurementEstimator *theEstimator;
  //const TransientTrackingRecHitBuilder *RHBuilder;
  //const KFTrajectorySmoother * theSmoother;
  // const KFTrajectoryFitter * theFitter;
  //  const TrackerGeometry * tracker;
  const MagneticField * magfield_;
  //TrajectoryStateTransform tsTransform;
  
  edm::ESWatcher<IdealMagneticFieldRecord> watchMagFieldRcd_;
  edm::ESWatcher<SiPixelTemplateDBObjectESProducerRcd> watchSiPixelTemplateRcd_;
  const SiPixelTemplateDBObject* templateDBobject_;
  std::vector<SiPixelTemplateStore> thePixelTemp_;

  edm::EDGetTokenT<TrajTrackAssociationCollection> t_trajTrack;
  
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomEsToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  edm::ESGetToken<SiPixelTemplateDBObject, SiPixelTemplateDBObjectESProducerRcd> siPixelTemplateEsToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoPerEventEsToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomPerEventEsToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
};

#endif
