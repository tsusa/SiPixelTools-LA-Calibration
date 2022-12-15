#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "CalibTracker/Records/interface/SiPixelTemplateDBObjectESProducerRcd.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelTemplateDBObject.h"
#include "CondFormats/SiPixelTransient/interface/SiPixelTemplateDefs.h"
#include "CondFormats/SiPixelTransient/interface/SiPixelTemplate.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyMap.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "../interface/SiPixelLorentzAngleAnalyser.h"
#include "../interface/FitFPixMuH.h"

using namespace std;
using namespace edm;
using namespace reco;

SiPixelLorentzAngleAnalyser::SiPixelLorentzAngleAnalyser(edm::ParameterSet const& conf) : 
  year_(conf.exists("year")?
        conf.getParameter<int>("year"): 2022),
  //  filename_(conf.exists("fileName")? 
  //	    conf.getParameter<std::string>("fileName"):""),
  filenameFit_(conf.getParameter<std::string>("fileNameFit")),
  source_(conf.exists("source")?
          conf.getParameter<std::string>("source"): "RECO"),
  subdetector_(conf.exists("subdetector")?
               conf.getParameter<std::string>("subdetector"): "BPix"),
  analysisType_(conf.exists("analysisType")?
                conf.getParameter<std::string>("analysisType"): "GA"),
  fitRegionMin_(conf.exists("fitRegionMin")?
                conf.getParameter<double>("fitRegionMin"): 5.),
  fitRegionMax_(conf.exists("fitRegionMax")?
                conf.getParameter<double>("fitRegionMax"): 280.),
  fitRegionLinMin_(conf.exists("fitRegionLinMin")?
                   conf.getParameter<double>("fitRegionLinMin"): 50.),
  fitRegionLinMax_(conf.exists("fitRegionLinMax")?
                   conf.getParameter<double>("fitRegionLinMax"): 235.),
  nBinsX_(conf.exists("nBinsX")? conf.getParameter<int>("nBinsX"): 200),
  xMin_(conf.exists("xMin")? conf.getParameter<double>("xMin"): -1000),
  xMax_(conf.exists("xMax")? conf.getParameter<double>("xMax"): 1000),
  nBinsY_(conf.exists("nBinsY")? conf.getParameter<int>("nBinsY"): 50),
  yMin_(conf.exists("yMin")? conf.getParameter<double>("yMin"): -100),
  yMax_(conf.exists("yMax")? conf.getParameter<double>("yMax"): 400),
  ptmin_(conf.exists("ptMin")? conf.getParameter<double>("ptMin"):3.),
  normChi2Max_(conf.exists("normChi2Max")?
               conf.getParameter<double>("normChi2Max"): 2.0),
  clustSizeYMin_(conf.exists("clustSizeYMin")?
                 conf.getParameter<vector<int> >("clustSizeYMin"):
                 std::vector({4,4,4,2})),
  clustSizeYMax_(conf.exists("clustSizeYMax")?
                 conf.getParameter<int>("clustSizeYMax"): 11),
  clustSizeXMax_(conf.exists("clustSizeXMax")?
                 conf.getParameter<int>("clustSizeXMax"): 5),
  residualMax_(conf.exists("residualMax")?
               conf.getParameter<double>("residualMax"): 0.005),
  clustChargeMax_(conf.exists("clustChargeMax")?
                  conf.getParameter<double>("clustChargeMax"): 50.),
  largePixCut_(conf.exists("largePixCut")?
               conf.getParameter<bool>("largePixCut"): true),
  maxNumberOfPixelsCut_(conf.exists("maxNumberOfPixelsCut")?
                        conf.getParameter<int>("maxNumberOfPixelsCut"): 25)
{  
  std::cout << "=========>>>>> Constructor" << std::endl;
  std::cout << "--------------------------------------"<< std::endl;
  std::cout << "FileNameFit: " << filenameFit_ << std::endl;
  std::cout << "Histogram Parameters: "<< std::endl;
  std::cout << "X bins: " << nBinsX_ << " "
            << xMin_ << " " << xMax_ << std::endl;
  std::cout << "Y bins: " << nBinsY_ << " "
            << yMin_ << " " << yMax_ << std::endl;
  std::cout << "--------------------------------------"<< std::endl;
  std::cout << "Fit Regions:"<< std::endl;
  std::cout << "fitRegionMin: " << fitRegionMin_ << std::endl;
  std::cout << "fitRegionMax: " << fitRegionMax_ << std::endl;
  std::cout << "fitRegionMinLin: " << fitRegionLinMin_ << std::endl;
  std::cout << "fitRegionMaxLin: " << fitRegionLinMax_ << std::endl;
  std::cout << "--------------------------------------"<< std::endl;
  std::cout << "Cut values: "<< std::endl;
  std::cout << "ptmin_ = " << ptmin_ << std::endl;
  std::cout << "normChi2Max_ = " << normChi2Max_ << std::endl;
  std::cout << "clustSizeYMin_ = ("
            << clustSizeYMin_[0] << ",  "
            << clustSizeYMin_[1] << ",  "
            << clustSizeYMin_[2] << ",  "
            << clustSizeYMin_[3] << ")" << std::endl;
  std::cout << "clustSizeYMax_ = " << clustSizeYMax_ << std::endl;
  std::cout << "clustSizeXMax_ = " << clustSizeXMax_ << std::endl;
  std::cout << "clustChargeMax_ = " << clustChargeMax_ << std::endl;
  std::cout << "residualMax_ = " << residualMax_ << std::endl;
  std::cout << "largePixCut_ = " << largePixCut_ << std::endl;
  std::cout << "maxNumberOfPixelsCut_ = "
            << maxNumberOfPixelsCut_ << std::endl;
  std::cout << "--------------------------------------"<< std::endl;
  usesResource("TFileService");
  
  if (source_  == "RECO"){
    t_trajTrack = consumes<TrajTrackAssociationCollection> (conf.getParameter<edm::InputTag>("src"));
    trackerHitAssociatorConfig_ = consumesCollector();
    siPixelTemplateEsToken_ = esConsumes<edm::Transition::BeginRun>();        
    geomPerEventEsToken_ = esConsumes();
  }
  
  geomEsToken_ = esConsumes<edm::Transition::BeginRun>();
  topoToken_ = esConsumes<edm::Transition::BeginRun>();
  topoPerEventEsToken_ = esConsumes();
  magFieldToken_ = esConsumes<MagneticField, IdealMagneticFieldRecord, 
 			      edm::Transition::BeginRun>();
  
  if (source_ == "Ntuple"){          
    std::string listFile(conf.getUntrackedParameter<std::string>("listFile", ""));
    std::cout << "==== >>>>>>>>> listFile: " << listFile << std::endl;
    if (listFile != ""){
      std::string line;
      ifstream ifp(listFile);
      std::cout << "is open: "<< ifp.is_open() << std::endl;
      if (ifp.is_open()){ 
        while (getline (ifp,line)){
          std::cout << "line: " << line << std::endl;
          if (line.length() < 10)
            continue;
          root_fileNames_.push_back(line);          
        }
        ifp.close();
      }
    }
  } // end if source == "Ntuple"  
} // end constructor

// Virtual destructor needed.
SiPixelLorentzAngleAnalyser::~SiPixelLorentzAngleAnalyser() {}  

void SiPixelLorentzAngleAnalyser::beginJob()
{  

  std::cout << "Start begin job" << std::endl;
  h_ptBC_= fs->make<TH1F>("h_ptBC", ";p_{T} (GeV);Entries", 200, 0, 50);
  h_ptAC_ = fs->make<TH1F>("h_ptAC", ";p_{T} (GeV);Entries", 200, 0, 50);
  h_ndofBC_ = fs->make<TH1F>("h_ndofBC", ";Ndof ;Entries", 51, -0.5, 50.5);
  h_ndofAC_ = fs->make<TH1F>("h_ndofAC", ";Ndof ;Entries", 51, -0.5, 50.5);
  h_chi2BC_ = fs->make<TH1F>("h_chi2BC", ";#Chi^{2} ;Entries", 100, 0, 100);
  h_chi2AC_ = fs->make<TH1F>("h_chi2AC", ";#Chi^{2} ;Entries", 100, 0, 100);
  h_normChi2BC_ = fs->make<TH1F>("h_normChi2BC", ";#Chi^{2}/Ndof ;Entries",
                           100, 0, 5);
  h_normChi2AC_ = fs->make<TH1F>("h_normChi2AC", ";#Chi^{2}/Ndof ;Entries",
				 100, 0, 5);
  h_nPixBC_ = fs->make<TH1F>("h_nPixBC", ";N_{pixels};Entries", 30, -0.5, 29.5);
  h_nPixAC_ = fs->make<TH1F>("h_nPixAC", ";N_{pixels};Entries", 30, -0.5, 29.5);
  h_cotBetaBC_ = fs->make<TH1F>("h_cotBetaBC",";Cot(#beta);Entries",
                          100, -10, 10);
  h_cotBetaAC_ = fs->make<TH1F>("h_cotBetaAC",";Cot(#beta);Entries",
                          100, -10, 10);
  h_cotAlphaBC_ = fs->make<TH1F>("h_cotAlphaBC",";Cot(#alpha);Entries",
                          100, -10, 10);
  h_cotAlphaAC_ = fs->make<TH1F>("h_cotAlphaAC",";Cot(#alpha);Entries",
                          100, -10, 10);


  for (int i=0; i< nLayersMax_; ++i){
    h_cotBetaLayerBC_[i] =
      fs->make<TH1F>(("h_cotBetaLayerBC_"+ std::to_string(i+1)).c_str(),
               ";Cot(#beta);Entries", 100, -10, 10);
    h_cotBetaLayerAC_[i] =
      fs->make<TH1F>(("h_cotBetaLayerAC_"+ std::to_string(i+1)).c_str(),
               ";Cot(#beta);Entries", 100, -10, 10);
  }
  h_largePixBC_ = fs->make<TH1F>("h_largePixBC", "Large Pixel", 2, -0.5, 1.5);
  h_largePixAC_ = fs->make<TH1F>("h_largePixAC","Large Pixel", 2, -0.5, 1.5);
  h_clusterSizeYBC_ = fs->make<TH1F>("h_clusterSizeYBC",
				     "Cluster Size Y [pixel]; Entries",
                               26, -0.5, 25.5);
  h_clusterSizeYAC_ = fs->make<TH1F>("h_clusterSizeYAC",
                               "Cluster Size Y [pixel]; Entries",
                               26, -0.5, 25.5);

  for (int i=0; i< nLayersMax_; ++i){
    h_clusterSizeYLayerBC_[i] =
      fs->make<TH1F>(("h_clusterSizeYLayerBC_" + std::to_string(i+1)).c_str(),
               ";Cluster Size Y [pixel]",  26, -0.5, 25.5);
    h_clusterSizeYLayerAC_[i] =
      fs->make<TH1F>(("h_clusterSizeYLayerAC_" + std::to_string(i+1)).c_str(),
               ";Cluster Size Y [pixel]",  26, -0.5, 25.5);
  }
  h_clusterSizeXBC_ = fs->make<TH1F>("h_clusterSizeXBC",
                               ";Cluster Size X [pixel]", 11, -0.5, 10.5);
  h_clusterSizeXAC_ = fs->make<TH1F>("h_clusterSizeXAC",
                               ";Cluster Size X [pixel]", 11, -0.5, 10.5);
  h_clusterChargeBC_ = fs->make<TH1F>("h_clusterChargeBC",
                                "Cluster charge (ke);Entries",
                                100, 0., 250.);
  h_clusterChargeAC_ = fs->make<TH1F>("h_clusterChargeAC",
				 "Cluster charge (ke);Entries",
				 100, 0., 250.);
  h_residualBC_ = fs->make<TH1F>("h_residualBC", ";residual (#cm);Entries",
			    200, 0, 0.1);
  h_residualAC_ = fs->make<TH1F>("h_residualAC",";residual (#cm);Entries",
				 200, 0, 0.1);
  
  //==========================================================================
  if (subdetector_ == "FPix"){         
    for (int i=0; i < nRingsMax_;++i){
      for (int j=0; j < nPanelsMax_;++j){
        for (int k=0; k < nZSidesMax_;++k){
          std::string name = "R"+std::to_string(i+1) +
            "_P" + std::to_string(j+1) + "_z" + std::to_string(k);
          std::string title = "Ring " + std::to_string (i+1) +
            "_Panel " + std::to_string(j+1) + "_z";
          title += (k == 0 ? ", z-": ", z+");

          mcshFPixAlpha_[i][j][k].Book(fs, name + "_test", "alpha", title, nBinsX_, nBinsY_);
          mcshFPixBeta_[i][j][k].Book(fs, name + "_test", "beta", title, nBinsX_, nBinsY_);
	  
          countMagField_[i][j][k] = 0;
          for (int m=0; m < 3;++m){
            std::string name_mgf = name + "_mfcomp_" +std::to_string(m);
            std::string title_mgf = title + "_mfcomp_" +std::to_string(m);

            magFieldHisto_[i][j][k][m] = fs->make<TH1D> (name_mgf.c_str(),
							 title_mgf.c_str(),
							 10000, -5., 5.);
            meanMagField_[i][j][k][m] = 0;
          }
        }
      }
    }
  } // if FPix and MCS
  
  //===========================================================================
  //mcshTest_.Book(fs, "test", "alpha", "test", 60, 60);  
  // just for some expaining plots
  // fixMe (some of the variables below should be removed)
  event_counter_ = 0;
  trackEventsCounter_ = 0;
  // 	trackcounter_ = 0;
  hitCounter_ = 0;
  usedHitCounter_ = 0;
  pixelTracksCounter_ = 0;
  
  nlay = 0;
  for (int i = 0; i < 4; i++){
     nModules_[i] = 0;
  }

  std::cout << "End begin job" << std::endl;
} // end beginJob

void SiPixelLorentzAngleAnalyser::beginRun(const edm::Run & iRun, const edm::EventSetup & es){
  
  std::cout << "Start begin run" << std::endl;
  if (source_ == "RECO"){
    const TrackerGeometry* geom = &es.getData(geomEsToken_);
    const TrackerTopology* tTopo = &es.getData(topoToken_);
    if (watchMagFieldRcd_.check(es)) {      
      magfield_ = &es.getData(magFieldToken_);
    }
    
    if (watchSiPixelTemplateRcd_.check(es)) {
      templateDBobject_ = &es.getData(siPixelTemplateEsToken_);
      if (!SiPixelTemplate::pushfile(*templateDBobject_, thePixelTemp_)) {
        edm::LogError("SiPixelLorentzAnglePCLWorker")
	  << "Templates not filled correctly. Check the sqlite file. Using SiPixelTemplateDBObject version "
          << (*templateDBobject_).version();
      }
    }    
  }
  
  /*
    PixelTopologyMap map = PixelTopologyMap(geom,tTopo);
    
    nlay = geom->numberOfLayers(PixelSubdetector::PixelBarrel);
    
    for (int i = 0; i < nlay; i++){
      nModules_[i] = map.getPXBModules(i + 1);
      }
    */
    // book histograms
  // char name[128];  
  std::cout << "End begin run" << std::endl;
}


// Functions that gets called by framework every event
void SiPixelLorentzAngleAnalyser::analyze(const edm::Event& e, const edm::EventSetup& es)
{    

  std::cout << "Start analyze" << std::endl;
  if (source_ == "Ntuple"){
    if (subdetector_ == "BPix"){
      std::cout << "Bpix Ntuple analysis" << std::endl;
      //ntupleAnalysisBpix();
      return;
    }
    if (subdetector_ == "FPix"){
      std::cout << "Calling nytupleAnalysisFPix" << std::endl;
      ntupleAnalysisFpix();
      std::cout << "Done fpix aanalysis...We should now return " << std::endl;
      return;
    }
  }
  
  event_counter_++;
  //std::cout << "Event counter: " << event_counter_ << std::endl;
  
  SiPixelTemplate templ(thePixelTemp_);  
  const TrackerTopology* const tTopo = &es.getData(topoPerEventEsToken_);
  const TrackerGeometry* tracker = &es.getData(geomPerEventEsToken_);
 

  TrackerHitAssociator* associate(0);
  
  // eout << associate << endl;
  // reset values
  module_=-1;
  layer_=-1;
  ladder_ = -1;
  isflipped_ = -1;
  pt_ = -999;
  eta_ = 999;
  phi_ = 999;
  pixinfo_.npix = 0;

  run_       = e.id().run();
  event_     = e.id().event();
  lumiblock_ = e.luminosityBlock();
  bx_        = e.bunchCrossing();
  orbit_     = e.orbitNumber();

  
  edm::Handle<TrajTrackAssociationCollection> trajTrackCollectionHandle;
  e.getByToken(t_trajTrack,trajTrackCollectionHandle);
  
  if(trajTrackCollectionHandle->size() >0){
    trackEventsCounter_++;
    if (trackEventsCounter_%1000 == 0)
      std::cout << "Looping over tracks: " << trackEventsCounter_ << std::endl;
    
    for(TrajTrackAssociationCollection::const_iterator it = trajTrackCollectionHandle->begin(); it!=trajTrackCollectionHandle->end();++it){
      const Track&      track = *it->val;
      const Trajectory& traj  = *it->key;
      
      // get the trajectory measurements
      std::vector<TrajectoryMeasurement> tmColl = traj.measurements(); 
      // 			TrajectoryStateOnSurface tsos = tsoscomb( itTraj->forwardPredictedState(), itTraj->backwardPredictedState() );
      pt_ = track.pt();
      eta_ = track.eta();
      phi_ = track.phi();
      chi2_ = traj.chiSquared();
      ndof_ = traj.ndof();
      //   if(pt_ < ptmin_) continue;
      // iterate over trajectory measurements
      std::vector<PSimHit> matched;
      //  h_tracks_->Fill(0);
      bool pixeltrack = false;
      for(std::vector<TrajectoryMeasurement>::const_iterator itTraj = tmColl.begin(); itTraj != tmColl.end(); itTraj++) {
	if(! itTraj->updatedState().isValid()) continue;
	
	TransientTrackingRecHit::ConstRecHitPointer recHit = itTraj->recHit();
	if(! recHit->isValid() || recHit->geographicalId().det() != DetId::Tracker ) continue;
	unsigned int subDetID = (recHit->geographicalId().subdetId());
	if( subDetID == PixelSubdetector::PixelBarrel || subDetID == PixelSubdetector::PixelEndcap){
	  
	  if(!pixeltrack){
	    //fixMe h_tracks_->Fill(1);
	    pixelTracksCounter_++;
	    std::cout << "pixelTrackCounter: " << pixelTracksCounter_ << std::endl;
	  }
	  pixeltrack = true;
	}
	if( subDetID == PixelSubdetector::PixelBarrel){
	  hitCounter_++;
				
	  DetId detIdObj = recHit->geographicalId();
	  const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (tracker->idToDet(detIdObj));
	  if(!theGeomDet) continue;
	  
	  const PixelTopology * topol = &(theGeomDet->specificTopology());
	  
	  if(!topol) continue;
	  
	  layer_ = tTopo->pxbLayer(detIdObj);
	  ladder_ = tTopo->pxbLadder(detIdObj);
	  module_ = tTopo->pxbModule(detIdObj);
	  float tmp1 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,0.)).perp();
	  float tmp2 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,1.)).perp();
	  if (tmp2<tmp1 ) isflipped_ = 1;
	  else isflipped_ = 0;
	  const SiPixelRecHit * recHitPix = dynamic_cast<const SiPixelRecHit *>((*recHit).hit());
	  if(!recHitPix) continue;
	  rechit_.x  = recHitPix->localPosition().x();
	  rechit_.y  = recHitPix->localPosition().y();
	  SiPixelRecHit::ClusterRef const& cluster = recHitPix->cluster();	
	
	  // fill entries in clust_
	  clust_.x = (cluster)->x();
	  clust_.y = (cluster)->y();
	  clust_.charge = (cluster->charge())/1000.;
	  clust_.size_x = cluster->sizeX();
	  clust_.size_y = cluster->sizeY();
	  clust_.maxPixelCol = cluster->maxPixelCol();
	  clust_.maxPixelRow = cluster->maxPixelRow();
	  clust_.minPixelCol = cluster->minPixelCol();
	  clust_.minPixelRow = cluster->minPixelRow();

	  // fill entries in pixinfo_:
	  fillPix(*cluster ,topol, pixinfo_);
	  // fill the trackhit info
	  TrajectoryStateOnSurface tsos=itTraj->updatedState();
	  if(!tsos.isValid()){
	    cout << "tsos not valid" << endl;
	    continue;	
	  }	
	  LocalVector trackdirection=tsos.localDirection();
	  LocalPoint trackposition=tsos.localPosition();
	
	  if(trackdirection.z()==0) continue;				
	  // the local position and direction
	  trackhit_.alpha = atan2(trackdirection.z(),trackdirection.x());
	  trackhit_.beta = atan2(trackdirection.z(),trackdirection.y());
	  trackhit_.gamma = atan2(trackdirection.x(),trackdirection.y());
	  trackhit_.x = trackposition.x();
	  trackhit_.y = trackposition.y();
	  
	  // get qScale_ = templ.qscale() and  templ.r_qMeas_qTrue();
	  
	  float cotalpha = 1./TMath::Tan(trackhit_.alpha);
          float cotbeta = 1./TMath::Tan(trackhit_.beta);
          float locBx = 1.;
          if(cotbeta < 0.)
            locBx = -1.;
          float locBz = locBx;
          if(cotalpha < 0.)
            locBz = -locBx;
	  
          auto  detId = detIdObj.rawId();
          int TemplID = -9999;
          TemplID = templateDBobject_->getTemplateID(detId);
          templ.interpolate(TemplID, cotalpha, cotbeta, locBz, locBx);
	  qScale_ = templ.qscale();
          rQmQt_ = templ.r_qMeas_qTrue();
	  
          // Surface deformation
          LocalPoint lp_track;        
          LocalPoint lp_rechit;
	  surface_derormation(topol, tsos, recHitPix, lp_track, lp_rechit); 

          rechitCorr_.x = lp_rechit.x();
          rechitCorr_.y = lp_rechit.y();
          trackhitCorrX_ = lp_track.x();
          trackhitCorrY_ = lp_track.y();
	  
	} // if BPix
	else if (subDetID == PixelSubdetector::PixelEndcap) {
	  DetId detIdObj = recHit->geographicalId();
	  const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> ( tracker->idToDet(detIdObj) );
	  if(!theGeomDet) continue;
	  const PixelTopology * topol = &(theGeomDet->specificTopology());
	  
	  if(!topol) continue;
	  sideF_ = tTopo->pxfSide(detIdObj);
	  diskF_ = tTopo->pxfDisk(detIdObj);
	  bladeF_ = tTopo->pxfBlade(detIdObj);
	  panelF_ = tTopo->pxfPanel(detIdObj);
	  moduleF_ = tTopo->pxfModule(detIdObj);
	  
	  const SiPixelRecHit * recHitPix = dynamic_cast<const SiPixelRecHit *>((*recHit).hit());
	  if(!recHitPix) continue;
	  rechitF_.x  = recHitPix->localPosition().x();
	  rechitF_.y  = recHitPix->localPosition().y();
	  SiPixelRecHit::ClusterRef const& cluster = recHitPix->cluster();	
	  // fill entries in clust_
	  clustF_.x = (cluster)->x();
	  clustF_.y = (cluster)->y();
	  clustF_.charge = (cluster->charge())/1000.;
	  clustF_.size_x = cluster->sizeX();
	  clustF_.size_y = cluster->sizeY();
	  clustF_.maxPixelCol = cluster->maxPixelCol();
	  clustF_.maxPixelRow = cluster->maxPixelRow();
	  clustF_.minPixelCol = cluster->minPixelCol();
	  clustF_.minPixelRow = cluster->minPixelRow();
	  // fill entries in pixinfo_:
	  fillPix(*cluster ,topol, pixinfoF_);
	  // fill the trackhit info
	  TrajectoryStateOnSurface tsos=itTraj->updatedState();
	  if(!tsos.isValid()){
	    cout << "tsos not valid" << endl;
	    continue;	
	  }	
	  LocalVector trackdirection=tsos.localDirection();	  
	  LocalPoint trackposition=tsos.localPosition();	 
	  if(trackdirection.z()==0) continue;				
	  // the local position and direction	
	  trackhitF_.alpha = atan2(trackdirection.z(),trackdirection.x());
	  trackhitF_.beta = atan2(trackdirection.z(),trackdirection.y());
	  trackhitF_.gamma = atan2(trackdirection.x(),trackdirection.y());
	  
	  trackhitF_.x = trackposition.x();
	  trackhitF_.y = trackposition.y();
	  
	  // get qScale_ = templ.qscale() 
	  float cotalpha = 1./TMath::Tan(trackhitF_.alpha);
          float cotbeta = 1./TMath::Tan(trackhitF_.beta);
          float locBx = 1.;
          if(cotbeta < 0.)
            locBx = -1.;
          float locBz = locBx;
          if(cotalpha < 0.)
            locBz = -locBx;
	  
	  LocalVector Bfield = theGeomDet->surface().toLocal(magfield_->inTesla(theGeomDet->surface().position()));
	  
	  magFieldF_[0] = Bfield.x();
	  magFieldF_[1] = Bfield.y();
	  magFieldF_[2] = Bfield.z();
	  
          auto  detId = detIdObj.rawId();
          int TemplID = -9999;
          TemplID = templateDBobject_->getTemplateID(detId);
          templ.interpolate(TemplID, cotalpha, cotbeta, locBz, locBx);
	  qScaleF_ = templ.qscale();
          rQmQtF_ = templ.r_qMeas_qTrue();
	  
          // Surface deformation
          LocalPoint lp_track;        
          LocalPoint lp_rechit;
	  surface_derormation(topol, tsos, recHitPix, lp_track, lp_rechit); 
	  rechitCorrF_.x = lp_rechit.x();
          rechitCorrF_.y = lp_rechit.y();
          trackhitCorrXF_ = lp_track.x();
          trackhitCorrYF_ = lp_track.y();
	  fpixMinimumClusterSize();
	}
      }	//end iteration over trajectory measurements
    } //end iteration over trajectories
  } // trajTrackCollection
} // end analyze

void SiPixelLorentzAngleAnalyser::endJob()  
{
  
  // FPix Minimal Cluster Size analysis  
  if (analysisType_ == "MCS" && subdetector_ == "FPix"){
    fstream my_file;
    my_file.open("results.txt", ios::out);
    
    for (int i=0; i < nRingsMax_;++i){
      for (int j=0; j < nPanelsMax_;++j){
	FitFPixMuH fitMuH;
	FPixMuH fmh;
        for (int k=0; k < nZSidesMax_;++k){
	  int retAlpha = mcshFPixAlpha_[i][j][k].Fit(fitRegionMax_);	  
	  int retBeta = mcshFPixBeta_[i][j][k].Fit(fitRegionMax_);
	  if (retAlpha || retBeta)
	    continue;
	  if (countMagField_[i][j][k] != 0){
	    for (int m=0; m < 3;++m){
              meanMagField_[i][j][k][m] /= countMagField_[i][j][k];
	    }
	    my_file << "\n";
	    fmh.bfield[0] = meanMagField_[i][j][k][0];
            fmh.bfield[1] = meanMagField_[i][j][k][1];
            fmh.bfield[2] = meanMagField_[i][j][k][2];
            fmh.shiftx =  mcshFPixAlpha_[i][j][k].GetMinCotAngle();
            fmh.shiftx_err =  mcshFPixAlpha_[i][j][k].GetMinCotAngleError();
            fmh.shifty =  mcshFPixBeta_[i][j][k].GetMinCotAngle();
            fmh.shifty_err =  mcshFPixBeta_[i][j][k].GetMinCotAngleError();
	    fitMuH.Add(fmh);	    
          }	  
        } // loop over z- sides
	
	fitMuH.Fit();
	std::cout << "Fitting Ring Panel: " <<  i+1 << " " << j+1 << " " <<  
	  fitMuH.GetMuH() << " " << fitMuH.GetMuHErr() << std::endl;
        std::cout << "---------------------------" << std::endl;
	my_file << "Fitting Ring Panel: " <<  i+1 << " " << j+1 <<  " " << 
	  fitMuH.GetMuH() << " " << fitMuH.GetMuHErr() << std::endl;
      } // loop over panels

    } // loop over Rings
    my_file.close();    
  } // if FPix       
} // endJobs

 inline void SiPixelLorentzAngleAnalyser::fillPix(const SiPixelCluster & LocPix, const PixelTopology * topol, Pixinfo& pixinfo)
   
{
  const std::vector<SiPixelCluster::Pixel>& pixvector = LocPix.pixels();
  pixinfo.npix = 0;
  for(std::vector<SiPixelCluster::Pixel>::const_iterator itPix = pixvector.begin(); itPix != pixvector.end(); itPix++){
    // 	for(pixinfo.npix = 0; pixinfo.npix < static_cast<int>(pixvector.size()); ++pixinfo.npix) {
    pixinfo.row[pixinfo.npix] = itPix->x;
    pixinfo.col[pixinfo.npix] = itPix->y;
    pixinfo.adc[pixinfo.npix] = itPix->adc;
    LocalPoint lp = topol->localPosition(MeasurementPoint(itPix->x + 0.5, itPix->y+0.5));
    pixinfo.x[pixinfo.npix] = lp.x();
    pixinfo.y[pixinfo.npix]= lp.y();
    pixinfo.npix++;
  }
}
/*
//fixMe
void SiPixelLorentzAngleAnalyser::findMean(int i, int i_ring)
{
double nentries = 0;
	
h_drift_depth_adc_slice_->Reset("ICE");
	
// determine sigma and sigma^2 of the adc counts and average adc counts
//loop over bins in drift width
for( int j = 1; j<= hist_drift_; j++){
if(_h_drift_depth_noadc_[i_ring]->GetBinContent(j, i) >= 1){
double adc_error2 = (_h_drift_depth_adc2_[i_ring]->GetBinContent(j,i) - _h_drift_depth_adc_[i_ring]->GetBinContent(j,i)*_h_drift_depth_adc_[i_ring]->GetBinContent(j, i) / _h_drift_depth_noadc_[i_ring]->GetBinContent(j, i)) /  _h_drift_depth_noadc_[i_ring]->GetBinContent(j, i);
_h_drift_depth_adc_[i_ring]->SetBinError(j, i, sqrt(adc_error2));
double error2 = adc_error2 / (_h_drift_depth_noadc_[i_ring]->GetBinContent(j,i) - 1.);
_h_drift_depth_[i_ring]->SetBinError(j,i,sqrt(error2));
} 
else{
_h_drift_depth_[i_ring]->SetBinError(j,i,0);
_h_drift_depth_adc_[i_ring]->SetBinError(j, i, 0);
}
      
}
*/
void SiPixelLorentzAngleAnalyser::surface_derormation(const PixelTopology *topol, 
						      TrajectoryStateOnSurface &tsos, 
						      const SiPixelRecHit *recHitPix,
						      LocalPoint &lp_track, 
						      LocalPoint &lp_rechit){

  LocalPoint trackposition=tsos.localPosition();	  
  const LocalTrajectoryParameters& ltp = tsos.localParameters();
  const Topology::LocalTrackAngles  localTrackAngles (ltp.dxdz(), ltp.dydz());
   
  std::pair<float, float> pixels_track = topol->pixel(trackposition, localTrackAngles);
  std::pair<float, float> pixels_rechit = topol->pixel(recHitPix->localPosition(),
						       localTrackAngles);
   
  lp_track = topol->localPosition(MeasurementPoint(pixels_track.first,
						   pixels_track.second));
   
  lp_rechit = topol->localPosition(MeasurementPoint(pixels_rechit.first,
						    pixels_rechit.second));
}
 
//===============================================================================
int SiPixelLorentzAngleAnalyser::fpixMinimumClusterSize(){
  
  //std::cout << "Running FPix Minimum ClusterSize" << std::endl;
  
  if (fabs(fabs(magFieldF_[2])-3.5) > 0.5)
    return 1;
     	 
  //--- pt cut
  h_ptBC_->Fill(pt_);
  if (pt_< ptmin_)
    return 1;
  h_ptAC_->Fill(pt_);
  
  //--- normChi2 cut 
  h_ndofBC_->Fill(ndof_);
  h_chi2BC_->Fill(chi2_);
  if (ndof_ == 0)
    return 1;  
  double chi2_ndof = chi2_/ndof_;
  h_normChi2BC_->Fill(chi2_ndof);
  if (chi2_ndof >= normChi2Max_)
    return 1;
  h_ndofAC_->Fill(ndof_);
  h_chi2AC_->Fill(chi2_);
  h_normChi2AC_->Fill(chi2_ndof);

  //--- large pixel cut
  bool large_pix = false;
  for (int j = 0; j <  pixinfoF_.npix; j++){
    int colpos = static_cast<int>(pixinfoF_.col[j]);
    if (pixinfoF_.row[j] == 0 || pixinfoF_.row[j] == 79 ||
        pixinfoF_.row[j] == 80 || pixinfoF_.row[j] == 159 ||
        colpos % 52 == 0 || colpos % 52 == 51 ){
      large_pix = true;
    }
  }

  h_largePixBC_->Fill(large_pix);
  if (largePixCut_ && large_pix)
    return 1;
  h_largePixAC_->Fill(large_pix);
  
  
  //--- cluster charge cut 
  
  h_clusterChargeBC_->Fill(clustF_.charge);
  if (clustF_.charge >= clustChargeMax_)      
    return 1;    
  h_clusterChargeAC_->Fill(clustF_.charge);
  
  //--- residual cut 
  double residual = sqrt(pow(trackhitCorrXF_ - rechitCorrF_.x, 2) +
			 pow(trackhitCorrYF_ - rechitCorrF_.y, 2));
    
  
  h_residualBC_->Fill(residual);
  if (residual > residualMax_)  
    return 1;
  h_residualAC_->Fill(residual);
  
  
  //------------------------------------
  
  int ringIdx = bladeF_ <= 22? 0: 1;
  int panelIdx = panelF_-1;
  int sideIdx = sideF_-1;
  double cotanAlpha  = TMath::Tan(TMath::Pi()/2. - trackhitF_.alpha);
  double cotanBeta  = TMath::Tan(TMath::Pi()/2. - trackhitF_.beta);
  
  countMagField_[ringIdx][panelIdx][sideIdx]++;
  for (int mfIdx=0;mfIdx<3;++mfIdx){
    magFieldHisto_[ringIdx][panelIdx][sideIdx][mfIdx]->Fill(magFieldF_[mfIdx]);
    meanMagField_[ringIdx][panelIdx][sideIdx][mfIdx] += magFieldF_[mfIdx];  }
  
  h_clusterSizeXBC_->Fill(clustF_.size_x);
  h_clusterSizeYBC_->Fill(clustF_.size_y);
  
  if (clustF_.size_y >= 2){
    h_clusterSizeYAC_->Fill(clustF_.size_y);
    mcshFPixAlpha_[ringIdx][panelIdx][sideIdx].Fill(cotanAlpha,
                                                        clustF_.size_x);
  }
  
  if (clust_.size_x  >= 0){
    h_clusterSizeXAC_->Fill(clustF_.size_x);
    mcshFPixBeta_[ringIdx][panelIdx][sideIdx].Fill(cotanBeta,
                                                       clustF_.size_y);
  }
  
  return 0;
}
//=====================================================================
void SiPixelLorentzAngleAnalyser::ntupleAnalysisFpix(){
  
  for (auto f_name: root_fileNames_){
    std::cout << f_name << std::endl;
    TFile f(f_name.c_str(),"READ");
    if(f.IsZombie()) {
      std::cout << "File " << f_name << " is zombie. "  <<  std::endl;
      f.Close();
      continue;      
    } // check if zombie
    
    //f.cd("lorentzAngle");
    
    std::cout << "Get Tree" << std::endl;
    TTree *T = (TTree*) f.Get("lorentzAngle/SiPixelLorentzAngleTreeForward_");
    std::cout <<  T << std::endl;
    
    T->SetBranchAddress("run", &run_);
    T->SetBranchAddress("event", &event_);
    T->SetBranchAddress("module", &moduleF_);
    T->SetBranchAddress("side", &sideF_);
    T->SetBranchAddress("blade", &bladeF_);
    T->SetBranchAddress("panel", &panelF_);
    //T->SetBranchAddress("isflipped", &isflipped_);
    T->SetBranchAddress("pt", &pt_);
    T->SetBranchAddress("eta", &eta_);
    T->SetBranchAddress("phi", &phi_);
    T->SetBranchAddress("chi2", &chi2_);
    T->SetBranchAddress("ndof", &ndof_);
    T->SetBranchAddress("trackhit", &trackhitF_);
    T->SetBranchAddress("npix", &pixinfoF_.npix);
    T->SetBranchAddress("rowpix", pixinfoF_.row);
    T->SetBranchAddress("colpix", pixinfoF_.col);
    T->SetBranchAddress("adc", pixinfoF_.adc);
    T->SetBranchAddress("xpix", pixinfoF_.x);
    T->SetBranchAddress("ypix", pixinfoF_.y);
    T->SetBranchAddress("clust", &clustF_); // charge is given in 1000 e
    T->SetBranchAddress("rechit", &rechitF_);
    T->SetBranchAddress("rechit_corr", &rechitCorrF_);
    T->SetBranchAddress("trackhitcorr_x", &trackhitCorrXF_);
    T->SetBranchAddress("trackhitcorr_y", &trackhitCorrYF_);
    T->SetBranchAddress("mf_x", &magFieldF_[0]);
    T->SetBranchAddress("mf_y", &magFieldF_[1]);
    T->SetBranchAddress("mf_z", &magFieldF_[2]);
    
    unsigned int n_entries = T->GetEntries();
    std::cout << "# of entries: " << n_entries << std::endl;
    for (unsigned int i_en = 0;  i_en < n_entries; ++i_en){
      if (i_en % 10000 == 0)
        std::cout << "Entry: " << i_en << std::endl;
      T->GetEntry(i_en);
      fpixMinimumClusterSize();
    } // loop over tree entries
  } // loop over root files
} // end ntupleAnalysis

//========================================================================


  
  
