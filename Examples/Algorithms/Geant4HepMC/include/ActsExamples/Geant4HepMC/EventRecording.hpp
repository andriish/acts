// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <mutex>
#include <G4VUserDetectorConstruction.hh>

class G4RunManager;

namespace ActsExamples {

class EventRecording final : public ActsExamples::BareAlgorithm {
 public:
  /// @class Config
  struct Config {
    std::string inputParticles = "";
    std::string outputHepMcTracks = "geant-outcome-tracks";

    std::unique_ptr<G4VUserDetectorConstruction> detectorConstruction = nullptr;

    /// random number seed 1
    int seed1 = 12345;
    /// random number seed 2
    int seed2 = 45678;

    /// List of processes that can be combined to a single vertex
    std::vector<std::string> processFilter;
    /// List of processes that should be recorded
    std::vector<std::string> eventSelectionProcess;
    /// List to veto events with certain processes
    std::vector<std::string> eventRejectionProcess;
  };

  /// Constructor
  EventRecording(Config&& cnf, Acts::Logging::Level level);
  ~EventRecording();

  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const final override;

 private:
  /// The config object
  Config m_cfg;
  /// G4 run manager
  std::unique_ptr<G4RunManager> m_runManager;

  // has to be mutable; algorithm interface enforces object constness
  mutable std::mutex m_runManagerLock;
};
}  // namespace ActsExamples