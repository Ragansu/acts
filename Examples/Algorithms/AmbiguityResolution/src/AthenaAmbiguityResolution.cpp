// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

ActsExamples::AthenaAmbiguityResolution::AthenaAmbiguityResolution(
    std::string name, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm(name, lvl) {}


ActsExamples::ConstTrackContainer
ActsExamples::AthenaAmbiguityResolution::prepareOutputTrack(
    const ActsExamples::ConstTrackContainer& tracks,
    std::vector<std::size_t>& goodTracks) const {
  std::shared_ptr<Acts::ConstVectorMultiTrajectory> trackStateContainer =
      tracks.trackStateContainerHolder();
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  trackContainer->reserve(goodTracks.size());
  // Temporary empty track state container: we don't change the original one,
  // but we need one for filtering
  auto tempTrackStateContainer =
      std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer solvedTracks{trackContainer, tempTrackStateContainer};
  solvedTracks.ensureDynamicColumns(tracks);

  for (auto&& iTrack : goodTracks) {
    auto destProxy = solvedTracks.getTrack(solvedTracks.addTrack());
    auto srcProxy = tracks.getTrack(iTrack);
    destProxy.copyFrom(srcProxy, false);
    destProxy.tipIndex() = srcProxy.tipIndex();
  }

  ConstTrackContainer outputTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      trackStateContainer};
  return outputTracks;
}

std::vector<int> ActsExamples::AthenaAmbiguityResolution::simpleScore(
 const ActsExamples::ConstTrackContainer& tracks, std::vector<std::map<std::size_t, Counter>>& counterMaps) const {

  std::vector<int> trackScore;
  int iTrack = 0;  
  // std::vector<std::size_t> measurements;
  // std::vector<std::vector<std::size_t>> measurementsPerTrack;

  // Loop over all the trajectories in the events
  for (auto track : tracks){
    auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(
      tracks.trackStateContainer(), track.tipIndex());
    int score = 100;
    auto counterMap = m_counterMap;
    
    if (track.chi2() > 0 && track.nDoF() > 0) {
      score+= std::log10(1.0-(track.chi2()/track.nDoF())); // place holder
    }


    // detector score is determined by the number of hits/hole/outliers * hit/hole/outlier score
    // here so instead of calculating nhits/nholes/noutliers per volume, 
    // we just loop over all volumes and add the score.

    for (long unsigned int i = 0; i < trajState.measurementVolume.size(); ++i){
      auto detector_it = m_volumeMap.find(trajState.measurementVolume[i]);
      if(detector_it != m_volumeMap.end()){
        auto detector = detector_it->second;
        score+=detector.hitsScoreWeight;
        counterMap[detector.detectorId].nhits++;
      }
    }

    for (long unsigned int i = 0; i < trajState.holeVolume.size(); ++i){
      auto detector_it = m_volumeMap.find(trajState.holeVolume[i]);
      if(detector_it != m_volumeMap.end()){
        auto detector = detector_it->second;
        score+=detector.holesScoreWeight;
        counterMap[detector.detectorId].nholes++;
      }
      else{
        ACTS_INFO("Detector not found");
        ACTS_INFO("Detector ID: " << trajState.holeVolume[i]);
      } 
    }
    for (long unsigned int i = 0; i < trajState.outlierVolume.size(); ++i){
      auto detector_it = m_volumeMap.find(trajState.outlierVolume[i]);
      if(detector_it != m_volumeMap.end()){
        auto detector = detector_it->second;
        score+=detector.outliersScoreWeight;
        counterMap[detector.detectorId].noutliers++;
      }
    }
      // TODO: add scored based on eta and phi

    trackScore.push_back(score);
    // ACTS_INFO("Track " << iTrack << " score: " << score);

    // measurementPerTrack.push_back(track.measurements();

    counterMaps.push_back(counterMap);
    iTrack++;


  } // end of loop over tracks
    
  return trackScore;
}

// place holder for goodTracks algorithm
std::vector<std::size_t> 
ActsExamples::AthenaAmbiguityResolution::solveAmbiguity(
    const ActsExamples::ConstTrackContainer& tracks ,std::vector<int> trackScore, std::vector<std::map<std::size_t, Counter>>& counterMaps) const {
  
  for(long unsigned int i=0; i<5; ++i){
    ACTS_INFO("Track " << i << " score: " << trackScore[i]);
  }
  std::vector<std::size_t> cleanTracks = getCleanedOutTracks(tracks, counterMaps);

  ACTS_INFO("Number of tracks: " << tracks.size());

  ACTS_INFO("Number of clean tracks: " << cleanTracks.size());

  std::vector<std::size_t> goodTracks;

  ACTS_INFO("Min score: " << m_minScore);


  for(long unsigned int i=0; i<cleanTracks.size(); ++i){
    ACTS_INFO("Track " << i << " score: " << trackScore[cleanTracks[i]]);
    if (trackScore[cleanTracks[i]] > m_minScore){
      goodTracks.push_back(cleanTracks[i]);
    }
  }
  ACTS_INFO("Number of good tracks: " << goodTracks.size());
  return goodTracks;
}


std::vector<std::size_t> ActsExamples::AthenaAmbiguityResolution::getCleanedOutTracks(
    const ActsExamples::ConstTrackContainer& tracks,   std::vector<std::map<std::size_t, Counter>>& counterMaps) const {
  std::vector<std::size_t> cleanTracks;
    enum TsosTypes {
    // A measurement not yet used in any other track
    UnusedHit   = 1,
    // A measurement shared with another track
    SharedHit   = 2,
    // A hit that needs to be removed from the track
    RejectedHit = 3,
    // an outlier, to be copied in case
    Outlier     = 4,
    // other TSOS types to be copied in case
    OtherTsos   = 5
  };
  // Loop over all detectors


  std::vector<int> tsosType = std::vector<int>(tracks.size(), OtherTsos);

  // std::vector<std::vector<std::size_t>> measurementsPerTrack;
  // boost::container::flat_map<std::size_t,boost::container::flat_set<std::size_t>> tracksPerMeasurement;
  // std::size_t numberOfTracks = tracks.size();



  // for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
  //   for (auto iMeasurement : measurementsPerTrack[iTrack]) {
  //     tracksPerMeasurement[iMeasurement].insert(iTrack);
  //   }
  // }

  // std::vector<std::size_t> sharedMeasurementsPerTrack = std::vector<std::size_t>(tracks.size(), 0);
  // for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
  //   for (auto iMeasurement : measurementsPerTrack[iTrack]) {
  //     if (tracksPerMeasurement[iMeasurement].size() > 1) {
  //       ++sharedMeasurementsPerTrack[iTrack];
  //     }
  //   }
  // }

  int iTrack = 0;
  for (const auto& track : tracks) {

    auto counterMap = counterMaps[iTrack];

    // init array
    tsosType[iTrack] = OtherTsos;


    // if we do not have a measurement, we should just mark it

    auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(tracks.trackStateContainer(), track.tipIndex());
    bool TrkCouldBeAccepted = true;
 
    for(long unsigned int i = 0; i< trajState.measurementVolume.size(); ++i){
      auto detector_it = m_volumeMap.find(trajState.measurementVolume[i]);
      if(detector_it == m_volumeMap.end()){
        continue;
      }
      auto detector = detector_it->second;
      // ACTS_INFO ("---> Found summary information");

      // ACTS_INFO ("---> Detector ID: " << detector.detectorId);
      // ACTS_INFO ("---> Number of hits: " << counterMap[detector.detectorId].nhits);
      // ACTS_INFO ("---> Number of holes: " << counterMap[detector.detectorId].nholes);
      // ACTS_INFO ("---> Number of outliers: " << counterMap[detector.detectorId].noutliers);


      if (counterMap[detector.detectorId].nhits < detector.minHits){
        TrkCouldBeAccepted = false;
      }

      else if (counterMap[detector.detectorId].nholes > detector.maxHoles){
        TrkCouldBeAccepted = false;
      }

      else if (counterMap[detector.detectorId].noutliers > detector.maxOutliers){
        TrkCouldBeAccepted = false;
      }
      else
      {
        TrkCouldBeAccepted = true;
      }
      

    }
    if (TrkCouldBeAccepted){
      cleanTracks.push_back(iTrack);
      ACTS_INFO("Track " << iTrack << " is clean");
    }
    iTrack++;

  }
  return cleanTracks;
}

// place holder for goodTracks algorithm
// OLD CODE IGNORE






    // Loop over all tracks
  
  // bool TrkCouldBeAccepted        = true;
  // // some counters used in the logic
  // int  numUnused         = 0;
  // int  numTRT_Unused     = 0;
  // int  numShared         = 0;
  // int  numWeightedShared = 0;
  // bool thishasblayer     = false;
  // bool hassharedblayer   = false;
  // bool hassharedpixel    = false;
  // bool firstisshared     = true; // logic is that it is set to false if we find a first hit unshared

  // // let's remember the last 2 ROTs on the track
  // const Trk::RIO_OnTrack* lastrot       = nullptr;
  // const Trk::RIO_OnTrack* lastbutonerot = nullptr;
  // int                     lastrotindex  = 0;


  // // get all TSOS the track
  // ConstTrackContainer tsos = track.trackStateContainer();

  
  // ACTS_DEBUG ("Study new Track "<< tracks<<"\t , it has "<<tsos->size()<<"\t track states");
  // ACTS_DEBUG ("trackId "<< trackId <<", subtrackId "<<subtrackId);

  // // is this a track from the pattern or a fitted track ?
  // bool ispatterntrack = (tracks->info().trackFitter()==Trk::TrackInfo::Unknown);
  // if (ispatterntrack) {
  //   ACTS_DEBUG ("==> this is a pattern track, outliers are good hits (reintegration) !");
  // } else {
  //   ACTS_DEBUG ("==> this is a refitted track, so we can use the chi2 ! ");
  // }

  // // some pre-processing of the summary information, if available, needed for special cuts
  // int nPixelDeadSensor = -1;
  // int nSCTDeadSensor   = -1;
  // int npixel    = 0;
  // int npixholes = 0;
  // int nsctholes = 0;
  // bool isTRT    = false;
  // auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(
  //     tracks.trackStateContainer(), track.tipIndex());
  
  // if (trajState.nStates > 0) {
  //   ACTS_VERBOSE ("---> Found summary information");
  //   nPixelDeadSensor = 0 //temporarily disabled 
  //   nSCTDeadSensor   = 0 //temporarily disabled
  //   npixel           = trajState.nMeasurements
  //   npixholes        = trajState.nHoles //temporarily assigned to nHoles
  //   nsctholes        = trajState.nHoles // temporarily assigned to nHoles
  // }
  // // set nDeadSensors to 0 in case trkSummary wasn't called with HoleSearch
  // // (i.e. number of deadSensors not available)
  // if (nPixelDeadSensor == -1) nPixelDeadSensor = 0;
  // if (nSCTDeadSensor   == -1) nSCTDeadSensor   = 0;
  // ACTS_VERBOSE ("---> Number of dead si sensors: " << nPixelDeadSensor + nSCTDeadSensor);

  // possible classification of TSOS

  // create an array of types for each TSOS
  

  // loop over all TSOS and classify them



    // ok, let try to get the ROT then
    // const Trk::RIO_OnTrack* rot = dynamic_cast <const Trk::RIO_OnTrack*> (meas);
    // if (!rot) {
    //   // could be a Pseudo-Measurement ?
    //   const Trk::PseudoMeasurementOnTrack* pseudo = dynamic_cast <const Trk::PseudoMeasurementOnTrack*> (meas);
    //   if (pseudo){
    //     ACTS_VERBOSE ("-> Copy pseudo measurement");
    //   } else {
    //     ACTS_WARNING ("-> Measurement is not a pseudo measurment, not yet supported, try to copy !");
    //   }
    //   tsosType[index] = OtherTsos;
    //   continue;
    // }

    //
    // we have a TSOS with a measurement, keep analysing it
    //

    // let's get some information about the measurement
    // const Identifier& id = rot->identify();
  //   bool isPixel         = m_detID->is_pixel(id);
  //   bool isBlayer        = isPixel ? m_detID->is_blayer(id) : false;
  //   // bool isoutlier       = (*iTsos)->type(Trk::TrackStateOnSurface::Outlier);

  //   // do we have some b-layer hit here (sorting of TSOS is inside out)
  //   if (isBlayer && (!isoutlier || ispatterntrack)) thishasblayer = true;  // we may reintegrate outliers from pattern

  //   // special cut to remove problematic combinations in the pixels
  //   if ( isPixel && !thishasblayer && npixholes>0           &&  // a pixel hit, no b-layer, we do have pixel holes
  //        ( ( npixel==1 && !isoutlier ) ||                       // one pixel on track, it is not an additional outlier
  //          ( ispatterntrack && npixel==0 && isoutlier) )   ) {  // pattern track, no pixels, but an outlier (we may reintegrate it)
  //     ACTS_VERBOSE ("-> Special case, problematic single pixel hit on track, reject it !");
  //     tsosType[index]    = RejectedHit;
  //     // mark track as bad !
  //     TrkCouldBeAccepted = false;
  //     continue;
  //   }

  //   // do we have an outlier (an not a pattern track) ?
  //   // if ( (isoutlier && !ispatterntrack) || !(m_detID->is_indet(id)) ) {
  //   if ( (isoutlier && !ispatterntrack) || !(m_detID->is_indet(id)) ) {
  //     ACTS_VERBOSE ("-> Prd is outlier on a fitter track (or not InDet), copy it over");
  //     tsosType[index] = Outlier;
  //     continue;
  //   }

  //   // let's check if this is a shared hit (even if it is an outlier on a pattern track) ?
  //   if (!prd_to_track_map.isUsed(*(rot->prepRawData()))) {
  //     if ( !isoutlier ) {
  //       ACTS_VERBOSE ("-> Prd is unused, copy it over");
  //     } else {
  //       ACTS_VERBOSE ("-> Prd is outlier on a pattern track and is unused, copy it over");
  //     }

  //     tsosType[index] = UnusedHit;
  //     // increase some counters
  //     if (isTRT) numTRT_Unused++;
  //     else numUnused++;
  //     // remember if first hit is shared, we need that later
  //     if (numShared == 0) firstisshared = false;
  //     // remember the last 2 ROTs
  //     lastbutonerot = lastrot;
  //     lastrot       = rot;
  //     lastrotindex  = index;

  //     continue;
  //   }

  //   //
  //   // ok, we have a shared hit
  //   //

  //   // do we have an outlier and a pattern track, but the hit is shared, so reject it (we reintegrate it otherwise)
  //   if ( isoutlier && ispatterntrack ) {
  //     ACTS_VERBOSE ("-> Shared Prd is outlier on a pattern track, we do not want to reintegrate it, so reject it ");
  //     tsosType[index]    = RejectedHit;
  //     // mark track as bad !
  //     TrkCouldBeAccepted = false;
  //     continue;
  //   }

  //   int numberOfTracksWithThisPrd = sharedMeasurementsPerTrack[index]

  //   // see if we try keeping it as a shared hit ?
  //   if ( numberOfTracksWithThisPrd < m_maxTracksPerPRD  && // we do not allow to share with to many tracks
  //        score > m_minScoreShareTracks                  && // score needs to be good
  //        (!isPixel || npixholes<=0)                    ) { // we do not share pixels if there are holes in pixels

  //     // special treatment of share split pixel clusters...
  //     if (m_doPixelClusterSplitting && isPixel) {
  //       // get pixel cluster
  //       const InDet::PixelCluster* clus = dynamic_cast <const InDet::PixelCluster*> (rot->prepRawData());

  //       if ( !clus ) {
  //         ACTS_WARNING ("-----> cast to Pixel cluster failed, should not happen !");
  //         // mark track as bad !
  //         TrkCouldBeAccepted = false; // we have to remove at least one PRD
  //         continue;
  //       } else {

  //         // split clusters are not allowed to be shared at all, unless
  //         const Trk::ClusterSplitProbabilityContainer::ProbabilityInfo &splitProb = splitProbContainer.splitProbability(clus);
  //         if ( splitProb.isSplit() )  {
  //           ACTS_VERBOSE ("-----> Pixel cluster is split, reject shared hit !!!");
  //           tsosType[index]    = RejectedHit;
  //           // mark track as bad !
  //           TrkCouldBeAccepted = false; // we have to remove at least one PRD
  //           continue;
  //         }

  //         // is cluster compatible with being a shared cluster ?
  //         // A.S.: also a hack for the max size: allows large clusters that are exluded from the splitter to be shared
  //         //       needs isExcluded() flag in the future
  //         if (splitProb.splitProbability1() < m_sharedProbCut && clus->rdoList().size() <= size_t(m_maxSplitSize) ) {
  //           ACTS_VERBOSE ("-----> Pixel cluster is not compatible with being shared (splitProb = "
  //                            << splitProb.splitProbability1() << ") , reject shared hit !!!");
  //           tsosType[index]    = RejectedHit;
  //           // mark track as bad !
  //           TrkCouldBeAccepted = false; // we have to remove at least one PRD
  //           continue;
  //         }
  //       }
  //     }

  //     ACTS_VERBOSE ("---> Shared hit, but good track, let's enter hit in the list and try to keep it !");
  //     tsosType[index]   = SharedHit;
  //     numShared++;                             // increase counter
  //     numWeightedShared += (isPixel ? 2 : 1);  // increase counter

  //     // some additional bookkeeping, we need this later
  //     if (isPixel) {
  //       // ME: bugfix, either shared blayer or shared pixel
  //       if (isBlayer) hassharedblayer = true;
  //       else  hassharedpixel = true;
  //     }
  //     // remember the last 2 ROTs
  //     lastbutonerot = lastrot;
  //     lastrot       = rot;
  //     lastrotindex  = index;

  //     continue;
  //   }

  //   // ok, we failed to keep the share hit, let's reject it.
  //   ACTS_VERBOSE ("---> Share hit and bad track, reject it !");
  //   tsosType[index]    = RejectedHit;
  //   // mark track as bad !
  //   TrkCouldBeAccepted = false; // we have to remove at least one PRD
  //   index++;
  // }

  // // total number of hits with dead modules
  // int totalSiHits = numUnused + nPixelDeadSensor + nSCTDeadSensor;

  // // special cut, do not allow the last hit to be to far away or after to many holes
  // if (ispatterntrack                                                                              && // pattern track and
  //     totalSiHits > m_minHits                                                                     && // we have enough hits on the track
  //     (lastrot && lastbutonerot)                                                                  && // has enough ROTs
  //     (lastrot->globalPosition()-lastbutonerot->globalPosition()).mag()>1000*CLHEP::mm) { // distance cut
  //   ACTS_DEBUG ("Special cut on distance, reject last hit on track !");
  //   tsosType[lastrotindex] = RejectedHit;
  //   numUnused--; // update counter
  //   // mark track as bad !
  //   TrkCouldBeAccepted     = false;
  // }

  // // special cut, do not allow the last hit to be after to many holes
  // if (ispatterntrack                                                                              && // pattern track and
  //     totalSiHits > m_minHits                                                                     && // we have enough hits on the track
  //     nsctholes>3 ) { // too many holes cut
  //   ACTS_DEBUG ("Special cut on too many holes, reject last hit on track !");
  //   tsosType[lastrotindex] = RejectedHit;
  //   numUnused--; // update counter
  //   // mark track as bad !
  //   TrkCouldBeAccepted     = false;
  // }



  // // get chi2/NDF, if track is fitted
  // if ( !ispatterntrack ) {
  //   double trackchi2 = 0.;
  //   if  (track.nDoF()>0 )
  //     trackchi2 = track.chi2 / track.nDoF();

  //   // if we have share hits and this is a bad track, we reject it
  //   if ( numShared > 0 && !ispatterntrack && trackchi2 > 3 ) {
  //     ACTS_DEBUG ("Shared hits, we have a bad chi2 track, mark it as bad !");
  //     // mark track as bad !
  //     TrkCouldBeAccepted = false;
  //   }
  // }



  //
  // now see what to do with the track
  //

  // best case, we like this track and it passes this complex if statement

//   if ( TrkCouldBeAccepted                                                    && // we did not mark the track as bad before
//        ( !hassharedblayer || npixholes<=1 )                                  && // if blayer, at most 1 pixel hole
//        !hassharedpixel                                                       && // no shared pixel hits
//        ( ( totalSiHits >= m_minHits      && numShared == 0 )                 || // no shared and enough hits OR
//          ( totalSiHits >= m_minNotShared && numWeightedShared <= m_maxShared && // shared hits and enough unique hits and shared hits with good quality
//            ( totalSiHits + std::min(numShared,m_maxShared) ) >= m_minHits && score > m_minScoreShareTracks ) ) ) {
//     ACTS_DEBUG ("=> Suggest to keep track with "<<numShared<<" shared hits !");
//     return std::make_tuple(static_cast<Trk::Track *>(nullptr),true); // keep input track

//     // ok, failed that one, can we recover the track ?
//   } else if ( numTRT_Unused >= nCutTRT           && // TRT track or no TRT at all (back tracking)
//               ( totalSiHits >= m_minHits         || // we have enough hits OR
//                 ( totalSiHits >= m_minNotShared  && // shared hits and enough unique hits and shared hits with good quality
//                   totalSiHits+std::min(numShared,m_maxShared) >= m_minHits && score > m_minScoreShareTracks ) ) ) {
//     // catch, if this is cosmics, accept the incoming track
//     if (m_cosmics) {
//       ACTS_DEBUG ("=> Cosmics, accept input track even with shared hits");
//       return std::make_tuple(static_cast<Trk::Track *>(nullptr),true); // keep input track;
//     }

//     ACTS_VERBOSE ("Trying to recover track, allow for some shared hits is possible.");

//     // new TSOS vector
//     std::vector<ConstTrackContainer> newTSOS;

//     // counter for the weighted number of added shared hits
//     int cntIns = 0;

//     // loop over all TSOS (and types) and copy the good ones overiTsos
//     auto iTsos    = tsos->begin();
//     auto iTsosEnd = tsos->end();

//     for (int index = 0 ; iTsos != iTsosEnd ; ++iTsos,++index ) {

//       // remove rejected hits
//       if (tsosType[index] == RejectedHit) {
//         ACTS_VERBOSE ("-> Dropping rejected hit");

//       } else if (tsosType[index] != SharedHit ) {
//         ACTS_VERBOSE ("-> Copy good TSOS");
//         newTSOS.push_back(*iTsos);

//       } else if (cntIns >= m_maxShared) {
//         ACTS_VERBOSE ("-> Too many share hits, dropping outer hit(s)");

//       } else {
//         ACTS_VERBOSE ("-> Try to recover a shared hit");

//         // get measurment from TSOS
//         const Trk::MeasurementBase* meas = (*iTsos)->measurementOnTrack();
//         // get ROT from this TSOS
//         const Trk::RIO_OnTrack*      rot = dynamic_cast <const Trk::RIO_OnTrack*> (meas);

//         if (!rot) {
//           ACTS_VERBOSE ("Cast to RIO_OnTrack failed, should never happen !");
//           continue;
//         }

//         // is it a pixel cluster ?
//         bool isPixel = m_detID->is_pixel(rot->identify());

//         // find out how many tracks use this hit already
//         Trk::PRDtoTrackMap::ConstPrepRawDataTrackMapRange range = prd_to_track_map.onTracks(*(rot->prepRawData()));
//         int                           numberOfTracksWithThisPrd = std::distance(range.first,range.second);
//         ACTS_VERBOSE ("---> number of tracks with this shared Prd: " << numberOfTracksWithThisPrd << " maxtracks: " << m_maxTracksPerPRD);

//         // check if this newly shared hit would exceed the shared hits limit of the already accepted track (**)
//         int  iShared        = 0;
//         int  othernpixel    = 0;
//         bool otherhasblayer = false;
//         if ( numberOfTracksWithThisPrd == 1 ) {
//           // @TODO is it ensured that the prds of all tracks have been registered already ?
//           std::vector< const Trk::PrepRawData* > prdsToCheck = m_assoTool->getPrdsOnTrack(prd_to_track_map,*(range.first->second));
//           for (const Trk::PrepRawData* prd : prdsToCheck) {
//             if (prd_to_track_map.isShared(*prd))
//               ++iShared;
//             if (m_detID->is_pixel(prd->identify())) {
//               othernpixel++;
//               if (m_detID->is_blayer(prd->identify())) otherhasblayer=true;
//             }
//           }
//         }

//         // now decide what to do, can we keep the shared hit
//         if ( numberOfTracksWithThisPrd < m_maxTracksPerPRD            && // number of tracks sharing hit
//              score > m_minScoreShareTracks                            && // score cut
//              (!isPixel || !hassharedblayer || npixholes <= 0)         && // shared b-layer only if no pixel holes
//              ( iShared < m_maxShared || (isPixel && !firstisshared) ) && // do not allow other accepted track to exceed the shared limit, if first pixel hit is shared
//              (!isPixel || thishasblayer == otherhasblayer )           && // only allow shared pixel if both have blayer or both not
//              (!isPixel || npixel >= othernpixel )                    ) { // only allow shared pixel if new track as at least as many pixel hits

//           ACTS_VERBOSE ("---> Accepted hit shared with " << numberOfTracksWithThisPrd << " tracks !");
//           newTSOS.push_back(*iTsos);
//           numUnused++; // update counter

//           // update shared hit counter
//           cntIns += (isPixel ? 2 : 1);

//         } else {
//           ACTS_VERBOSE ("---> Reject hit shared with " << numberOfTracksWithThisPrd << " tracks !");
//         }
//       }
//     }

//     // this still may happen per (**) see above.
//     if ( numUnused+nPixelDeadSensor+nSCTDeadSensor < m_minHits || newTSOS.size() <= 3 ) {
//       ACTS_DEBUG ("=> Too few hits, reject track with shared hits");
//       return std::make_tuple(static_cast<Trk::Track *>(nullptr),false); // reject input track;
//     }

//     // check that this is not the input track
//     if ( newTSOS.size() == tsos->size() ) {
//       ACTS_DEBUG ("=> Recovered input track, accept it !");
//       return std::make_tuple(static_cast<Trk::Track *>(nullptr),true); // keep input track;
//     } else {
//       // ok, done, create subtrack
//       // Trk::Track* newTrack = createSubTrack(newTSOS,track.nDoF);
//       // if (!newTrack) {
//       //   ACTS_DEBUG ("=> Failed to create subtrack");
//       //   return std::make_tuple(static_cast<Trk::Track *>(nullptr),false); // reject input track;
//       // }

//       Trk::TrackInfo info;
//       info.addPatternRecoAndProperties(ptrTrack->info());
//       Trk::TrackInfo newInfo;
//       newInfo.setPatternRecognitionInfo(Trk::TrackInfo::InDetAmbiTrackSelectionTool);
//       info.addPatternReco(newInfo);
//       newTrack->info().addPatternReco(ptrTrack->info());

//       ACTS_DEBUG ("=> Successfully created subtrack with shared hits recovered !");
//       return std::make_tuple(newTrack,false); // new cleaned track, reject input track ;
//     }
//   } else {
//     ACTS_DEBUG ("=> Track is recommended to be dropped");
//   }

//   return std::make_tuple(static_cast<Trk::Track *>(nullptr),false); // reject input track;

// }
