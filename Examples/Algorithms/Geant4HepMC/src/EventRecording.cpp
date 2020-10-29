// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4HepMC/EventRecording.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include <iostream>
#include <stdexcept>
#include <FTFP_BERT.hh>
#include "EventAction.hpp"
#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hpp"
#include "RunAction.hpp"
#include "SteppingAction.hpp"

#include <HepMC3/GenParticle.h>

ActsExamples::EventRecording::~EventRecording() {
  m_runManager = nullptr;
}

ActsExamples::EventRecording::EventRecording(
    ActsExamples::EventRecording::Config&& cnf, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("EventRecording", level),
      m_cfg(std::move(cnf)),
      m_runManager(std::make_unique<G4RunManager>()) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }
  if (m_cfg.outputHepMcTracks.empty()) {
    throw std::invalid_argument("Missing output event collection");
  }
  if (!m_cfg.detectorConstruction) {
    throw std::invalid_argument("Missing detector construction object");
  }

  /// Now set up the Geant4 simulation
  m_runManager->SetUserInitialization(m_cfg.detectorConstruction.release());
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new ActsExamples::RunAction());
  m_runManager->SetUserAction(
      new ActsExamples::EventAction(m_cfg.processesCombine));
  m_runManager->SetUserAction(
      new ActsExamples::PrimaryGeneratorAction(m_cfg.seed1, m_cfg.seed2));
  m_runManager->SetUserAction(
      new ActsExamples::SteppingAction(m_cfg.processesReject));
  m_runManager->Initialize();
}

std::pair<double, double>
  ActsExamples::EventRecording::collectMaterialBackwards(HepMC3::ConstGenVertexPtr vertex) const
  {
	double x0 = 0.;
	double l0 = 0.;
	  
	HepMC3::ConstGenParticlePtr currentParticle = vertex->particles_in()[0];
	const int id = currentParticle->attribute<HepMC3::IntAttribute>("TrackID")->value();
	
	while(vertex && !vertex->particles_in().empty() && vertex->particles_in()[0]->attribute<HepMC3::IntAttribute>("TrackID")
	 && vertex->particles_in()[0]->attribute<HepMC3::IntAttribute>("TrackID")->value() == id)
	{
		currentParticle = vertex->particles_in()[0];
		const double stepLength = currentParticle->attribute<HepMC3::DoubleAttribute>("StepLength")->value();
		x0 += stepLength / currentParticle->attribute<HepMC3::DoubleAttribute>("NextX0")->value();
		l0 += stepLength / currentParticle->attribute<HepMC3::DoubleAttribute>("NextL0")->value();
		vertex = currentParticle->production_vertex();
	}
	return std::make_pair(x0, l0);
  }
  
std::pair<double, double>
ActsExamples::EventRecording::collectMaterialForwards(HepMC3::ConstGenVertexPtr vertex) const
{
	double x0 = 0.;
	double l0 = 0.;
	
	HepMC3::ConstGenParticlePtr currentParticle = vertex->particles_out()[0];
	const int id = currentParticle->attribute<HepMC3::IntAttribute>("TrackID")->value();
	
	while(vertex && !vertex->particles_out().empty())
	{
		currentParticle = nullptr;
		for(const auto& particleOut : vertex->particles_out())
			if(particleOut->attribute<HepMC3::IntAttribute>("TrackID") && particleOut->attribute<HepMC3::IntAttribute>("TrackID")->value() == id)
			{
				currentParticle = particleOut;
				break;
			}
		const double stepLength = currentParticle->attribute<HepMC3::DoubleAttribute>("StepLength")->value();
		x0 += stepLength / currentParticle->attribute<HepMC3::DoubleAttribute>("NextX0")->value();
		l0 += stepLength / currentParticle->attribute<HepMC3::DoubleAttribute>("NextL0")->value();
		vertex = currentParticle->end_vertex();
	}
	
	return std::make_pair(x0, l0);
}
  
ActsExamples::ProcessCode ActsExamples::EventRecording::execute(
    const ActsExamples::AlgorithmContext& context) const {
  // ensure exclusive access to the geant run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);

  // Retrieve the initial particles
  const auto initialParticles =
      context.eventStore.get<ActsExamples::SimParticleContainer>(
          m_cfg.inputParticles);

  // Storage of events that will be produced
  std::vector<HepMC3::GenEvent> events;
  events.reserve(initialParticles.size());

  // Storage of interaction distances
  std::vector<std::tuple<bool, double, double>> interactionDistances;
  
  for (const auto& part : initialParticles) {
    // Prepare the particle gun
    ActsExamples::PrimaryGeneratorAction::instance()->prepareParticleGun(part);

    // Begin with the simulation
    m_runManager->BeamOn(1);

    // Test if the event was aborted
    if (SteppingAction::instance()->eventAborted()) {
      continue;
    }

	// Set event start time
    HepMC3::GenEvent event = ActsExamples::EventAction::instance()->event();
    HepMC3::FourVector shift(0., 0., 0., part.time() / Acts::UnitConstants::mm);
    event.shift_position_by(shift);
    
    // Set beam particle properties
    HepMC3::FourVector beamMom4(part.momentum4()[0], part.momentum4()[1], part.momentum4()[2], part.momentum4()[3]);
    auto beamParticle = event.particles()[0];
    beamParticle->set_momentum(beamMom4);
    beamParticle->set_pid(part.pdg());

    // Set beam particle properties
    const Acts::Vector4D momentum4 =
        part.momentum4() / Acts::UnitConstants::GeV;
    HepMC3::FourVector beamMom4(momentum4[0], momentum4[1], momentum4[2],
                                momentum4[3]);
    auto beamParticle = event.particles()[0];
    beamParticle->set_momentum(beamMom4);
    beamParticle->set_pid(part.pdg());

    if (m_cfg.processSelect.empty()) {
      // Store the result
      events.push_back(std::move(event));
    } else {
      bool storeEvent = false;
      // Test if the event has a process of interest in it
      for (const auto& vertex : event.vertices()) {
        if (vertex->id() == -1) {
          vertex->add_particle_in(beamParticle);
        }
        const std::vector<std::string> vertexAttributes =
            vertex->attribute_names();
        for (const auto& att : vertexAttributes) {
          if ((vertex->attribute_as_string(att).find(m_cfg.processSelect) !=
               std::string::npos) &&
              !vertex->particles_in().empty() &&
              vertex->particles_in()[0]->attribute<HepMC3::IntAttribute>(
                  "TrackID") &&
              vertex->particles_in()[0]
                      ->attribute<HepMC3::IntAttribute>("TrackID")
                      ->value() == 1) {
            storeEvent = true;
            break;
          }
        }
        if (storeEvent) {
          break;
        }
      }
      // Store the result
      if (storeEvent) {
        // Remove vertices without outgoing particles
        for (auto it = event.vertices().crbegin();
             it != event.vertices().crend(); it++) {
          if ((*it)->particles_out().empty()) {
            event.remove_vertex(*it);
          }
        }
        events.push_back(std::move(event));
      }
      // Store the result
      if(storeEvent)
      {
		  			  // Remove vertices without outgoing particles
			  for(auto it = event.vertices().crbegin(); it != event.vertices().crend(); it++)
			  {
				if((*it)->particles_out().empty())
				{
					event.remove_vertex(*it);
				}
				}
		events.push_back(std::move(event));
		
	}
	else
	{
		if(m_cfg.recordInteractionDistances)
		{
			auto material = storeEvent ? collectMaterialBackwards(processVertex) : collectMaterialForwards(event.vertices()[0]);
			interactionDistances.push_back(std::make_tuple(storeEvent, material.first, material.second));
		}
	}
    }
  }

  ACTS_INFO(initialParticles.size() << " initial particles provided");
  ACTS_INFO(events.size() << " tracks generated");

  // Write the recorded material to the event store
  context.eventStore.add(m_cfg.outputHepMcTracks, std::move(events));
  if(m_cfg.recordInteractionDistances)
	context.eventStore.add(m_cfg.outputInteraction, std::move(interactionDistances));

  return ActsExamples::ProcessCode::SUCCESS;
}
