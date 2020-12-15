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
#include <HepMC3/GenParticle.h>

#include "EventAction.hpp"
#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hpp"
#include "RunAction.hpp"
#include "SteppingAction.hpp"

#include <HepMC3/GenParticle.h>
#include <Eigen/Eigenvalues> 
#include <random>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include <TH1F.h>
#include <TCanvas.h>
#include "Acts/Definitions/Units.hpp"
#include <TFile.h>

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
  
  if(m_cfg.covarianceSample)
  {    
	  Acts::BoundVector stddev = Acts::BoundVector::Zero();
	  stddev[Acts::eBoundLoc0] = 15. * Acts::UnitConstants::um;
	  stddev[Acts::eBoundLoc1] = 80. * Acts::UnitConstants::um;
	  stddev[Acts::eBoundTime] = 0.1 * Acts::UnitConstants::ns;
	  stddev[Acts::eBoundPhi] = 1. * Acts::UnitConstants::degree;
	  stddev[Acts::eBoundTheta] = 1.5 * Acts::UnitConstants::degree;
	  stddev[Acts::eBoundQOverP] = 1. * Acts::UnitConstants::e / 10. * Acts::UnitConstants::GeV;
	  
	  Acts::BoundSymMatrix corr = Acts::BoundSymMatrix::Identity();
	  //~ corr(Acts::eBoundLoc0, Acts::eBoundLoc1) = corr(Acts::eBoundLoc1, Acts::eBoundLoc0) = 0.125;
	  //~ corr(Acts::eBoundLoc0, Acts::eBoundPhi) = corr(Acts::eBoundPhi, Acts::eBoundLoc0) = 0.25;
	  //~ corr(Acts::eBoundLoc1, Acts::eBoundTheta) = corr(Acts::eBoundTheta, Acts::eBoundLoc1) = -0.25;
	  //~ corr(Acts::eBoundTime, Acts::eBoundQOverP) = corr(Acts::eBoundQOverP, Acts::eBoundTime) = 0.125;
	  //~ corr(Acts::eBoundPhi, Acts::eBoundTheta) = corr(Acts::eBoundTheta, Acts::eBoundPhi) = -0.25;
	  //~ corr(Acts::eBoundPhi, Acts::eBoundQOverP) = corr(Acts::eBoundPhi, Acts::eBoundQOverP) = -0.125;
	  //~ corr(Acts::eBoundTheta, Acts::eBoundQOverP) = corr(Acts::eBoundTheta, Acts::eBoundQOverP) = 0.5;

	  m_cfg.covariance = stddev.asDiagonal() * corr * stddev.asDiagonal();

	Eigen::SelfAdjointEigenSolver<Acts::BoundSymMatrix> eigenSolver(m_cfg.covariance);

	m_transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }
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

  /// Build the histograms       
  std::vector<TH1F*> histos;
  histos.resize(Acts::eBoundSize);
  for(unsigned int i = 0; i < Acts::eBoundSize; i++)
  {
	  histos[i] = new TH1F("", "", 100, -5, 5);
  }
  ///

  for (const auto& part : initialParticles) {
	  for(unsigned int i = 0; i < m_cfg.numSamples; i++)
	  {
		ACTS_VERBOSE("Starting variation number " << i);
		 // Prepare the particle gun
		 if(m_cfg.covarianceSample)
		  {
			  ActsExamples::PrimaryGeneratorAction::instance()->prepareParticleGun(sampleFromCovariance(part, histos));
		  } else {
			ActsExamples::PrimaryGeneratorAction::instance()->prepareParticleGun(part);
		}
		
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
		const Acts::Vector4D momentum4 =
			part.fourMomentum() / Acts::UnitConstants::GeV;
		HepMC3::FourVector beamMom4(momentum4[0], momentum4[1], momentum4[2],
									momentum4[3]);
		auto beamParticle = event.particles()[0];
		beamParticle->set_momentum(beamMom4);
		beamParticle->set_pid(part.pdg());

		if (m_cfg.processSelect.empty()) {
			// Add the beam particle
		  for (const auto& vertex : event.vertices()) {
//~ std:cout << "Vtx ID: " << vertex->id() << std::endl;
				if (vertex->id() == -1) {
				  vertex->add_particle_in(beamParticle);
//~ std::cout << "Found" << std::endl;
				  break;
				}
			}
			
			// Remove vertices without outgoing particles
			for (auto it = event.vertices().crbegin();
				 it != event.vertices().crend(); it++) {
			  if ((*it)->particles_out().empty()) {
				event.remove_vertex(*it);
			  }
			}
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
		}
	}
  }
  
  TFile* tf = new TFile("G4InitialParameters.root", "RECREATE");
  TCanvas* tc = new TCanvas();
  for(unsigned int i = 0; i < Acts::eBoundSize; i++)
  {
	histos[i]->Draw();
	tc->Print(("Parameter_" + std::to_string(i) + ".png").c_str());
	histos[i]->Write(("InitialParameter" + std::to_string(i)).c_str());
	delete(histos[i]);
  }
  tf->Write();
  tf->Close();

  ACTS_INFO(initialParticles.size() << " initial particles provided");
  ACTS_INFO(events.size() << " tracks generated");

  // Write the recorded material to the event store
  context.eventStore.add(m_cfg.outputHepMcTracks, std::move(events));

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::SimParticle
ActsExamples::EventRecording::sampleFromCovariance(const ActsExamples::SimParticle& particle, const std::vector<TH1F*>& histos) const {
	  SimParticle result = particle;
	  
	  Acts::CurvilinearTrackParameters mean(particle.fourPosition(), particle.unitDirection(), particle.charge(), particle.absoluteMomentum());
	  
        static std::mt19937 gen{ std::random_device{}() };
        static std::normal_distribution<> dist;

        const Acts::BoundVector sample = mean.parameters() 
			+ m_transform * Acts::BoundVector{mean.parameters().size()}.unaryExpr(
				[&](ActsExamples::SimParticle::Scalar /*x*/) { return dist(gen); });
	
	for(unsigned int i = 0; i < Acts::eBoundSize; i++)
		histos[i]->Fill((sample[i] - mean.parameters()[i]) / std::sqrt(m_cfg.covariance(i, i))); // TODO: Pull only valid for diagonal matrix
	
	Acts::GeometryContext gctx;
	 const auto freeSample = Acts::detail::transformBoundToFreeParameters(mean.referenceSurface(), gctx, sample);

	result.setPosition4(freeSample.template segment<4>(Acts::eFreePos0));
	result.setDirection(freeSample.template segment<3>(Acts::eFreeDir0));
	result.setAbsoluteMomentum(std::abs(1. / freeSample[Acts::eFreeQOverP]));
	return result;
}