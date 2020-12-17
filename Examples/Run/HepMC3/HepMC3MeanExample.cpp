// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/HepMC/MeanCalculator.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Options.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <boost/program_options.hpp>

void addCustomOptions(::boost::program_options::options_description& desc) {
  //~ using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("mean-max-step-size", value<double>()->default_value(std::numeric_limits<double>::max()),
      "Maximum step size in propagation");
  opt("mean-max-steps", value<unsigned int>()->default_value(1000),
      "Maximum number of steps in propagation");
  opt("mean-tolerance", value<double>()->default_value(1e-4),
      "Maximum integration error tolerance in propagation");
}
 
ActsExamples::MeanCalculator::Config
readCustomOptions(const ::boost::program_options::variables_map& vars) {
	ActsExamples::MeanCalculator::Config cfg;
  cfg.maxStepSize = vars["mean-max-step-size"].as<double>();
  cfg.maxSteps = vars["mean-max-steps"].as<unsigned int>();
  cfg.tolerance = vars["mean-tolerance"].as<double>();
  return cfg;
}
      
int main(int argc, char** argv) {
	
	auto detector = std::make_shared<DD4hepDetector>();
	
  // Declare the supported program options.
  // Setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  addCustomOptions(desc);
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addInputOptions(desc);
  ActsExamples::Options::addHepMC3ReaderOptions(desc);
  ActsExamples::Options::addGeometryOptions(desc);
  detector->addOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

// TODO: Bfield? 

  // The detector
  auto [trackingGeometry, contextDecorators] = ActsExamples::Geometry::build(vm, *detector);

  // Create the readers
  // Read particles (initial states) and clusters from CSV files
  auto particleReader = ActsExamples::Options::readCsvParticleReaderConfig(vm);
  particleReader.inputStem = "particles";
  particleReader.outputParticles = "particles";

  // Read the HepMC events
  auto hepMC3ReaderConfig = ActsExamples::Options::readHepMC3ReaderOptions(vm);
  hepMC3ReaderConfig.outputEvents = "hepmc-events";

  // Create the algorithm      
  ActsExamples::MeanCalculator::Config meanConfig = readCustomOptions(vm);
  meanConfig.inputEvents = hepMC3ReaderConfig.outputEvents;
  meanConfig.inputParticles = particleReader.outputParticles;
  meanConfig.trackingGeometry = trackingGeometry;

  // Add to the sequencer
  sequencer.addReader(std::make_shared<ActsExamples::CsvParticleReader>(
      particleReader, logLevel));
  sequencer.addReader(std::make_shared<ActsExamples::HepMC3AsciiReader>(
      hepMC3ReaderConfig, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::MeanCalculator>(
      std::move(meanConfig), logLevel));

  // Run
  return sequencer.run();
}
