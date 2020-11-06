// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Io/Root/detail/NuclearInteractionParametrisation.hpp"

#include <mutex>
#include <vector>

namespace ActsExamples {

/// @brief This class takes fractions of recorded events that represent the
/// effect of a nuclear interaction and produces histograms and parameters which
/// can be used for a parametrisation based simultion of nuclear interaction.
/// Since the parameters are based on the set of all provided events, during the
/// event loop newly provided events are stored until the end of the run. Then
/// all parts are calculated and written to file.
class RootNuclearInteractionParametersWriter final
    : public WriterT<std::vector<
          std::tuple<ActsExamples::SimParticle, ActsExamples::SimParticle,
                     std::vector<ActsExamples::SimParticle>>>> {
 public:
  struct Config {
    /// Input collection to map measured hits to simulated hits.
    std::string inputEventFractions;
    /// output directory.
    std::string outputDir;
    /// output filename.
    std::string outputFilename = "parameters.root";
    /// file access mode.
    std::string fileMode = "RECREATE";

    /// Number of bins used for the interaction probability distributions
    unsigned int interactionProbabilityBins = 1e6;
    /// Number of bins used for the momentum distributions
    unsigned int momentumBins = 1e6;
    /// Number of bins used for the invariant mass distributions
    unsigned int invariantMassBins = 1e6;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  RootNuclearInteractionParametersWriter(const Config& cfg,
                                         Acts::Logging::Level lvl);
  ~RootNuclearInteractionParametersWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] event Fraction of an event that will be stored in @p
  /// m_eventFractionCollection
  ProcessCode writeT(
      const AlgorithmContext& /*ctx*/,
      const std::vector<
          std::tuple<ActsExamples::SimParticle, ActsExamples::SimParticle,
                     std::vector<ActsExamples::SimParticle>>>& event)
      final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  std::vector<NuclearInteractionParametrisation::EventFraction>
      m_eventFractionCollection;  ///< The recorded fractions of events
};
}  // namespace ActsExamples
