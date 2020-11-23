// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/Fatras/FatrasAlgorithm.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"
#include "ActsFatras/Kernel/Process.hpp"
#include "ActsFatras/Physics/StandardPhysicsLists.hpp"
#include "ActsExamples/Fatras/RootNuclearInteractionParametersReader.hpp"
#include "ActsFatras/Physics/NuclearInteraction/NuclearInteraction.hpp"
#include "ActsFatras/Physics/NuclearInteraction/Parameters.hpp"

#include <utility>

#include <boost/program_options.hpp>


class TH1F;



namespace ActsExamples {
namespace Options {

/// Add Fatras options.
///
/// @param desc The options description to add options to
void addNuclearInteractionOptions(Description& desc);

/// Reads the parametrisation and provides the parametrisation
ActsFatras::detail::MultiParticleParametrisation readParametrisations(
    const std::vector<std::string>& fileNames, const std::vector<int>& nSimulatedEvents);
      
/// Read Fatras options to create the algorithm config.
///
/// @tparam simulator_t type of the simulation kernel
/// @param vars         the variables to read from
/// @param simulator    the simulation kernel
template <typename simulator_t>
void readNuclearInteractionConfig(
    const Variables& variables, simulator_t& simulator) {

  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("NuclearInteractionOptions", Acts::Logging::INFO))

  
  const auto nuclearInteractionParametrisations = variables["fatras-nuclear-interaction-parametrisation"].as<std::vector<std::string>>();
  const auto nSimulatedEvents = variables["fatras-simulated-events-nuclear-interaction-parametrisation"].as<std::vector<int>>(); 
   	
	if(nuclearInteractionParametrisations.empty() && nSimulatedEvents.empty()) {
		ACTS_WARNING("No parametrisation for the nuclear interaction provided.");
		return;
	} else { 
		if(nuclearInteractionParametrisations.size() != nSimulatedEvents.size())
		{
			ACTS_WARNING("Unequal number of files and simulated provided: " 
					<< nuclearInteractionParametrisations.size() << " files and " << nSimulatedEvents.size() << " simulated events");
			return;
		}
	}	
	ACTS_VERBOSE(nuclearInteractionParametrisations.size() << " parametrisation files provided.");
	
	auto& chargedNuclearInteraction = simulator.charged.physics.template get<ActsFatras::detail::ParametrisedNuclearInteraction>();
	auto& neutralNuclearInteraction = simulator.neutral.physics.template get<ActsFatras::detail::ParametrisedNuclearInteraction>();
	
	const auto mpp = readParametrisations(nuclearInteractionParametrisations, nSimulatedEvents);
	ACTS_VERBOSE("Parametrisations for nuclear interaction from " << mpp.size() << " particles provided");
	
	chargedNuclearInteraction.physics = mpp;
	neutralNuclearInteraction.physics = mpp;
}

}  // namespace Options
}  // namespace ActsExamples