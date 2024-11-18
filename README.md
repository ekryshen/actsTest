# actsTest
Testing ACTS with simple telescope tracker

How to fit with modified mass hypothesis:
* Core/include/Acts/EventData/ParticleHypothesis.hpp
  inherits from
  Core/include/Acts/EventData/GenericParticleHypothesis.hpp

* Examples/Algorithms/TrackFitting/src/RefittingAlgorithm.cpp <- example on how to get input fitted tracks, fit and get output fitted tracks
  Important input: std::shared_ptr<TrackFitterFunction> fit;
  auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options, calibrator, surfSequence, tracks);

  const Acts::BoundTrackParameters initialParams(track.referenceSurface().getSharedPtr(), track.parameters(), track.covariance(), track.particleHypothesis());
  
  for (const auto& track : inputTracks) { // loop over tracks
    for (auto state : track.trackStatesReversed()) { // loop over track states

    }
  }


