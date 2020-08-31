// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ReconTruthBase.hpp"

#include <boost/program_options.hpp>

#include "ACTFW/Digitization/HitSmearing.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Fitting/FittingAlgorithm.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Csv/CsvOptionsReader.hpp"
#include "ACTFW/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ACTFW/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ACTFW/Io/Root/RootEffectiveTrajectoryWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/TruthTracking/ParticleSmearing.hpp"
#include "ACTFW/TruthTracking/TruthTrackFinder.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

using namespace Acts::UnitLiterals;

/// @brief Setup for the Truth Reconstruction
///
/// @param variables The boost variable map to resolve
/// @param sequencer The framework sequencer
/// @param tGeometry The TrackingGeometry for the tracking setup
/// truth reco
void
FW::setupReconTruth(
    const FW::Options::Variables&                 variables,
    FW::Sequencer&                                sequencer,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry)
{

  // Read the log level
  Acts::Logging::Level logLevel = FW::Options::readLogLevel(variables);
  // Output path
  auto outputDir
      = FW::ensureWritableDirectory(variables["output-dir"].as<std::string>());
  auto rnd = std::make_shared<FW::RandomNumbers>(
      FW::Options::readRandomNumbersConfig(variables));

  auto magneticField = FW::Options::readBField(variables);

  // Create smeared measurements
  FW::HitSmearing::Config hitSmearingCfg;
  hitSmearingCfg.inputSimulatedHits    = "hits";
  hitSmearingCfg.outputSourceLinks     = "sourcelinks";
  hitSmearingCfg.outputHitParticlesMap = "hit_particles_map";
  hitSmearingCfg.sigmaLoc0             = 25_um;
  hitSmearingCfg.sigmaLoc1             = 100_um;
  hitSmearingCfg.sigmaGlob0             = 200_um;
  hitSmearingCfg.sigmaGlob1             = 200_um;
  hitSmearingCfg.sigmaGlob2             = 200_um;
  hitSmearingCfg.randomNumbers         = rnd;
  hitSmearingCfg.trackingGeometry      = trackingGeometry;
  sequencer.addAlgorithm(
      std::make_shared<FW::HitSmearing>(hitSmearingCfg, logLevel));

  // The fitter needs the measurements (proto tracks) and initial
  // track states (proto states). The elements in both collections
  // must match and must be created from the same input particles.
  const auto& inputParticles = "particles_initial";
  // Create truth tracks
  FW::TruthTrackFinder::Config trackFinderCfg;
  trackFinderCfg.inputParticles       = inputParticles;
  trackFinderCfg.inputHitParticlesMap = hitSmearingCfg.outputHitParticlesMap;
  trackFinderCfg.outputProtoTracks    = "prototracks";
  sequencer.addAlgorithm(
      std::make_shared<FW::TruthTrackFinder>(trackFinderCfg, logLevel));
  // Create smeared particles states
  FW::ParticleSmearing::Config particleSmearingCfg;
  particleSmearingCfg.inputParticles        = inputParticles;
  particleSmearingCfg.outputTrackParameters = "smearedparameters";
  particleSmearingCfg.randomNumbers         = rnd;
  // Gaussian sigmas to smear particle parameters
  particleSmearingCfg.sigmaD0    = 20_um;
  particleSmearingCfg.sigmaD0PtA = 30_um;
  particleSmearingCfg.sigmaD0PtB = 0.3 / 1_GeV;
  particleSmearingCfg.sigmaZ0    = 20_um;
  particleSmearingCfg.sigmaZ0PtA = 30_um;
  particleSmearingCfg.sigmaZ0PtB = 0.3 / 1_GeV;
  particleSmearingCfg.sigmaPhi   = 1_degree;
  particleSmearingCfg.sigmaTheta = 1_degree;
  particleSmearingCfg.sigmaPRel  = 0.01;
  particleSmearingCfg.sigmaT0    = 1_ns;
  sequencer.addAlgorithm(
      std::make_shared<FW::ParticleSmearing>(particleSmearingCfg, logLevel));

  // setup the fitter
  FW::FittingAlgorithm::Config fitter;
  fitter.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  fitter.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  fitter.inputInitialTrackParameters
      = particleSmearingCfg.outputTrackParameters;
  fitter.outputTrajectories = "trajectories";
  fitter.fit                = FW::FittingAlgorithm::makeFitterFunction(
      trackingGeometry, magneticField, logLevel);
  sequencer.addAlgorithm(
      std::make_shared<FW::FittingAlgorithm>(fitter, logLevel));

  // write tracks from fitting
  FW::RootEffectiveTrajectoryWriter::Config trackWriter;
  trackWriter.inputParticles    = inputParticles;
  trackWriter.inputTrajectories = fitter.outputTrajectories;
  trackWriter.outputDir         = outputDir;
  trackWriter.outputFilename    = "tracks.root";
  trackWriter.outputTreename    = "tracks";
  sequencer.addWriter(
      std::make_shared<FW::RootEffectiveTrajectoryWriter>(trackWriter, logLevel));

  // write reconstruction performance data
  FW::TrackFinderPerformanceWriter::Config perfFinder;
  perfFinder.inputParticles       = inputParticles;
  perfFinder.inputHitParticlesMap = hitSmearingCfg.outputHitParticlesMap;
  perfFinder.inputProtoTracks     = trackFinderCfg.outputProtoTracks;
  perfFinder.outputDir            = outputDir;
  sequencer.addWriter(
      std::make_shared<FW::TrackFinderPerformanceWriter>(perfFinder, logLevel));
  FW::TrackFitterPerformanceWriter::Config perfFitter;
  perfFitter.inputParticles    = inputParticles;
  perfFitter.inputTrajectories = fitter.outputTrajectories;
  perfFitter.outputDir         = outputDir;
  sequencer.addWriter(
      std::make_shared<FW::TrackFitterPerformanceWriter>(perfFitter, logLevel));
}