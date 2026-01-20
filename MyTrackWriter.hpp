#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <string>
#include <vector>

using ActsExamples::AlgorithmContext;
using ActsExamples::ConstTrackContainer;
using ActsExamples::WriterT;
using ActsExamples::ProcessCode;

class TFile;
class TTree;

class MyTrackWriter final : public ActsExamples::WriterT<ConstTrackContainer> {
 public:
  struct Config {
    std::string inputTracks;
    std::string filePath = "tracks.root";
    std::string treeName = "tracks";
    std::string fileMode = "RECREATE";
  };

  MyTrackWriter(const Config& config, Acts::Logging::Level level);
  ~MyTrackWriter() override;
  ProcessCode finalize() override;
  const Config& config() const { return m_cfg; }

 protected:
  ProcessCode writeT(const AlgorithmContext& ctx, const ConstTrackContainer& tracks) override;

 private:
  Config m_cfg;
  TFile* m_outputFile{nullptr};
  TTree* m_outputTree{nullptr};
  std::uint32_t m_eventNr{0};
  std::vector<std::uint32_t> m_trackNr;
  std::vector<float> m_chi2Sum;
  std::vector<unsigned int> m_NDF;
  std::vector<float> m_eLOC0_fit;
  std::vector<float> m_eLOC1_fit;
  std::vector<float> m_ePHI_fit;
  std::vector<float> m_eTHETA_fit;
  std::vector<float> m_eQOP_fit;
  std::vector<std::vector<int>> m_measurementIds;
};
