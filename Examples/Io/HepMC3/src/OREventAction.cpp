// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "OREventAction.hpp"
#include <stdexcept>
#include "ORSteppingAction.hpp"
#include <G4Event.hh>
#include <G4RunManager.hh>

ActsExamples::OREventAction* ActsExamples::OREventAction::s_instance = nullptr;

ActsExamples::OREventAction*
ActsExamples::OREventAction::instance()
{
  // Static acces function via G4RunManager
  return s_instance;
}

ActsExamples::OREventAction::OREventAction() : G4UserEventAction()
{
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    s_instance = this;
  }
}

ActsExamples::OREventAction::~OREventAction()
{
  s_instance = nullptr;
}

void
ActsExamples::OREventAction::BeginOfEventAction(const G4Event*)
{
  ORSteppingAction::instance()->clear();
  m_event = std::make_shared<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);
}

void
ActsExamples::OREventAction::EndOfEventAction(const G4Event*)
{
	std::cout << "Particles: " << m_event->particles().size() << " | " << "Vertices: " << m_event->vertices().size() << std::endl;
	std::cout << "Number of steps: " << ORSteppingAction::instance()->counter() << std::endl;
	ORSteppingAction::instance()->counter() = 0;
	
	for(const auto& part : m_event->particles())
	{	
		if(!part->production_vertex())
			std::cout << "Production vertex missing" << std::endl;
		if(!part->end_vertex())
			std::cout << "End vertex missing" << std::endl;
	}
	for(const auto& vert : m_event->vertices())
	{
		auto proc = vert->attribute<HepMC3::StringAttribute>("Process");
		if(proc && (proc->value() == "Transportation" || proc->value() == "Death"))
			continue;
		if(vert->particles_out().size() > 1)
		std::cout << (proc ? proc->value() : "") << ": " << vert->particles_in().size() << "(" 
			<< (vert->particles_in().size() > 0 ? vert->particles_in()[0]->pid() : 0) << ") -> " << vert->particles_out().size() << std::endl;
	}
	std::cout << "Checks done" << std::endl;
}

void
ActsExamples::OREventAction::clear()
{
	m_event = nullptr;
	ORSteppingAction::instance()->clear();
}