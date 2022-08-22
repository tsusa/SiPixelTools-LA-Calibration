#ifndef CalibTracker_SiPixelLorentzAngle_SiPixelLorentzAngle_h
#define CalibTracker_SiPixelLorentzAngle_SiPixelLorentzAngle_h

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
#include "TROOT.h"


/**
* 
* Class to determine the lorentz angle in the barrel pixel detector 
* returns a txt file with the fit for every of the 8 rings in the 3 layers
* 
*/

// ggiurgiu@fnal.gov : remove namespace 12/27/09
//namespace
//{
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
//}

class SiPixelLorentzAngle : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> 
{
public:

  explicit SiPixelLorentzAngle(const edm::ParameterSet& conf);

  virtual ~SiPixelLorentzAngle();
  virtual void beginJob();
  virtual void endJob();
  void beginRun(const edm::Run& r, const edm::EventSetup& c) override;
  void endRun(const edm::Run& r, const edm::EventSetup& c) override {}  
  void analyze(const edm::Event& e, const edm::EventSetup& c) override;

 private:

  edm::Service<TFileService> fs;
  void fillPix(const SiPixelCluster & LocPix, const PixelTopology * topol, Pixinfo& pixinfo);
  void surface_derormation(const PixelTopology *topol, 
			   TrajectoryStateOnSurface &tsos, 
			   const SiPixelRecHit *recHitPix,
			   LocalPoint &lp_track, 
			   LocalPoint &lp_rechit);


  void findMean(int i, int i_ring);
  
  TTree* SiPixelLorentzAngleTree_;
  TTree* SiPixelLorentzAngleTreeForward_;
  
  // tree branches barrel
  int run_;
  long int event_;
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
  int eventF_;  
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
  
  // parameters from config file

  std::string filename_;
  std::string filenameFit_;
  double ptmin_;
  bool simData_;
  double normChi2Max_;
  int clustSizeYMin_;
  double residualMax_;
  double clustChargeMax_;
  int hist_depth_;
  int hist_drift_;

  TrackerHitAssociator::Config trackerHitAssociatorConfig_;
  // histogram etc
  int hist_x_;
  int hist_y_;
  double min_x_;
  double max_x_;
  double min_y_;
  double max_y_;
  double width_;
  double min_depth_;
  double max_depth_;
  double min_drift_;
  double max_drift_;
  
  std::map<int, TH2F*> _h_drift_depth_adc_;
  std::map<int, TH2F*> _h_drift_depth_adc2_;
  std::map<int, TH2F*> _h_drift_depth_noadc_;
  std::map<int, TH2F*> _h_drift_depth_;
  TH1F* h_drift_depth_adc_slice_;
  std::map<int, TH1F*> _h_mean_;
  TH2F *h_cluster_shape_adc_;
  TH2F *h_cluster_shape_noadc_;
  TH2F *h_cluster_shape_;
  TH2F *h_cluster_shape_adc_rot_;
  TH2F *h_cluster_shape_noadc_rot_;
  TH2F *h_cluster_shape_rot_;
  TH1F *h_tracks_;
  
  
  int event_counter_, trackEventsCounter_,pixelTracksCounter_, hitCounter_, usedHitCounter_;
  
  int nlay;
  int nModules_[4];
  
  // CMSSW classes needed
  PropagatorWithMaterial  *thePropagator;
  PropagatorWithMaterial  *thePropagatorOp;
  KFUpdator *theUpdator;
  Chi2MeasurementEstimator *theEstimator;
  const TransientTrackingRecHitBuilder *RHBuilder;
  const KFTrajectorySmoother * theSmoother;
  const KFTrajectoryFitter * theFitter;
  const TrackerGeometry * tracker;
  const MagneticField * magfield;
  TrajectoryStateTransform tsTransform;
 
  edm::ESWatcher<SiPixelTemplateDBObjectESProducerRcd> watchSiPixelTemplateRcd_;
  const SiPixelTemplateDBObject* templateDBobject_;
  std::vector<SiPixelTemplateStore> thePixelTemp_;

  edm::EDGetTokenT<TrajTrackAssociationCollection> t_trajTrack;
  edm::EDGetTokenT<edm::SimTrackContainer>  tok_simTk_;  


  edm::EDGetTokenT<edm::SimVertexContainer> tok_simVtx_;
  edm::EDGetTokenT<edm::PCaloHitContainer>  tok_caloEB_;
  edm::EDGetTokenT<edm::PCaloHitContainer>  tok_caloEE_;
  edm::EDGetTokenT<edm::PCaloHitContainer>  tok_caloHH_;
  
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomEsToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  edm::ESGetToken<SiPixelTemplateDBObject, SiPixelTemplateDBObjectESProducerRcd> siPixelTemplateEsToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoPerEventEsToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomPerEventEsToken_;
};

#endif
