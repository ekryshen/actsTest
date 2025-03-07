// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "MySpacePointMaker.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderOptions.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
//#include "ActsExamples/Utilities/GroupBy.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <utility>

ActsExamples::MySpacePointMaker::MySpacePointMaker(Config cfg, Acts::Logging::Level lvl)
  : IAlgorithm("MySpacePointMaker", lvl), m_cfg(std::move(cfg)) {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);

  auto spConstructor =
      [](const Acts::Vector3& pos, std::optional<double> t,
         const Acts::Vector2& cov, std::optional<double> varT,
         boost::container::static_vector<Acts::SourceLink, 2> slinks)
      -> SimSpacePoint {
    return SimSpacePoint(pos, t, cov[0], cov[1], varT, std::move(slinks));
  };

  auto spBuilderConfig = Acts::SpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = m_cfg.trackingGeometry;
  m_slSurfaceAccessor.emplace(IndexSourceLink::SurfaceAccessor{*m_cfg.trackingGeometry});
  spBuilderConfig.slSurfaceAccessor.connect<&IndexSourceLink::SurfaceAccessor::operator()>(&m_slSurfaceAccessor.value());
  m_spacePointBuilder = Acts::SpacePointBuilder<SimSpacePoint>(spBuilderConfig, spConstructor, Acts::getDefaultLogger("SpacePointBuilder", lvl));
}

ActsExamples::ProcessCode ActsExamples::MySpacePointMaker::execute(const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);

  // function to access measurement parameters using source links
  auto accessor = [&measurements](Acts::SourceLink slink) {
    const auto islink = slink.get<IndexSourceLink>();
    const ConstVariableBoundMeasurementProxy meas = measurements.getMeasurement(islink.index());
    return std::make_pair(meas.fullParameters(), meas.fullCovariance());
  };

  // fill front and back strip vectors of pairs(source link, pair<stripEnd1, stripEnd2>)
  std::vector<std::pair<Acts::SourceLink, std::pair<Acts::Vector3, Acts::Vector3>>> frontStrips;
  std::vector<std::pair<Acts::SourceLink, std::pair<Acts::Vector3, Acts::Vector3>>> backStrips;

  double rMax = 1300; // FIXME Extract from surface dimensions
  for (auto& isl : measurements.orderedIndices()) {
    const auto geoId = isl.geometryId();
    const auto volumeId = geoId.volume();
    const auto layerId = geoId.layer();
    ACTS_DEBUG("volumeId= " << volumeId << " layerId=" << layerId);
    Acts::SourceLink slink{isl}; 
    const auto [par, cov] = accessor(slink);
    const Acts::Surface* surface = m_slSurfaceAccessor.value()(slink);
    // TODO: more realistic strip dimensions including inner radii and half-station splitting
    auto gpos1 = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0],-rMax), Acts::Vector3());
    auto gpos2 = surface->localToGlobal(ctx.geoContext, Acts::Vector2(par[0], rMax), Acts::Vector3());
    if (layerId%4==2) {
      frontStrips.emplace_back(std::make_pair(slink, std::make_pair(gpos1, gpos2)));
    } else if (layerId%4==0) {
      backStrips.emplace_back(std::make_pair(slink, std::make_pair(gpos1, gpos2)));
    }
  }
  ACTS_DEBUG("making strip pairs:" << "back: " << frontStrips.size() << " front: " << backStrips.size());

  // make space points
  SimSpacePointContainer spacePoints;
  for (auto& fstrip : frontStrips) {
    float fz = fstrip.second.second[2];
    for (auto& bstrip : backStrips) {
      float bz = bstrip.second.second[2];
      if (fabs(fz - bz) > 20) continue;
      std::vector<Acts::SourceLink> slinks = {fstrip.first, bstrip.first};
      auto strippair = std::make_pair(fstrip.second, bstrip.second);
      Acts::SpacePointBuilderOptions spOpt{strippair, accessor};
      m_spacePointBuilder.buildSpacePoint(ctx.geoContext, slinks, spOpt, std::back_inserter(spacePoints));
    }
  }

  spacePoints.shrink_to_fit();

  ACTS_DEBUG("Created " << spacePoints.size() << " space points");
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ActsExamples::ProcessCode::SUCCESS;
}
