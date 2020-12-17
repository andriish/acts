// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/HepMC/Plot.hpp"

struct TrackSummary;

namespace ActsExamples {

Plot plot;

class MeanCalculator final : public ActsExamples::BareAlgorithm {
 public:
  /// @class Config
  struct Config {
    /// The input collection
    std::string inputEvents;
    /// The initial particles
    std::string inputParticles;
    
    /// The output collection
    std::string outputSimulationProcesses = "event-fraction";

    /// The detector
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    
    double maxStepSize = std::numeric_limits<double>::max();
    unsigned int maxSteps = 1000;
    double tolerance = 1e-4;
  };

  /// Constructor
  MeanCalculator(Config&& cnf, Acts::Logging::Level level);
  ~MeanCalculator();

  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const final override;

 private:
   
	void
	plotMean(const std::vector<TrackSummary>& summaries) const;
 
  /// The config object
  Config m_cfg;
  
};
}  // namespace ActsExamples
