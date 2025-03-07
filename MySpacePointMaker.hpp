#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {
struct AlgorithmContext;

class MySpacePointMaker final : public IAlgorithm {
 public:
  struct Config {
    std::string inputMeasurements;
    std::string outputSpacePoints;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::vector<Acts::GeometryIdentifier> geometrySelection;
  };
  MySpacePointMaker(Config cfg, Acts::Logging::Level lvl);
  ProcessCode execute(const AlgorithmContext& ctx) const override;
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  std::optional<IndexSourceLink::SurfaceAccessor> m_slSurfaceAccessor;
  Acts::SpacePointBuilder<SimSpacePoint> m_spacePointBuilder;
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};
  WriteDataHandle<SimSpacePointContainer> m_outputSpacePoints{this, "OutputSpacePoints"};
};
}
