// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Physics/Decay/Decay.hpp"

//~ // CLHEP
//~ #include "CLHEP/Units/SystemOfUnits.h"
//~ #include "CLHEP/Units/PhysicalConstants.h"
//~ #include "AtlasHepMC/GenParticle.h"


#include "G4RunManager.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "Acts/Utilities/Definitions.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"
#include "QGSP_BERT.hh"
#include "ActsFatras/Physics/Decay/G4DetectorConstruction.hpp"

/**AlgTool constructor for ParticleDecayHelper*/
ActsFatras::Decay::Decay() : m_g4RunManager(initG4RunManager()) {}
      
/** decay */
std::vector<ActsFatras::Particle> 
ActsFatras::Decay::decayParticle(const ActsFatras::Particle& parent) const {
  // return vector for children
  std::vector<Particle> children;

  int pdgCode = parent.pdg();
    
  G4ParticleDefinition* pDef = m_pdgToG4Conv.getParticleDefinition(pdgCode);
  if(!pDef)
  {
    return children;
  }
  
  G4DecayTable* dt = pDef->GetDecayTable();
  if(!dt)
  {
    return children;
  }
     
  G4VDecayChannel* channel = dt->SelectADecayChannel();
  if(!channel)
  {
    return children;
  }
  
  G4DecayProducts* products = channel->DecayIt();
  if(!products)
  {
    return children;
  }

  const Particle::Vector4 mom4 = parent.momentum4();
  products->Boost( mom4[Acts::eMom0] / mom4[Acts::eEnergy],
                   mom4[Acts::eMom1] / mom4[Acts::eEnergy],
                   mom4[Acts::eMom2] / mom4[Acts::eEnergy]);

  G4int nProducts = products->entries();
  for(G4int i = 0; i < nProducts; i++)
  {
    G4DynamicParticle* prod = products->PopProducts();
    if(!prod)
    {
      continue;
    }
    
    // the decay product
    const G4ThreeVector& mom= prod->GetMomentum();
    constexpr double convertEnergy = Acts::UnitConstants::GeV / CLHEP::GeV;
    Acts::Vector3D amgMom( mom.x(), mom.y(), mom.z() );
	amgMom *= convertEnergy;

	int32_t pdg = prod->GetPDGcode();
	Particle childParticle(Barcode(), static_cast<Acts::PdgParticle>(pdg));
	childParticle.setPosition4(parent.position4()).setAbsMomentum(amgMom.norm()).setDirection(amgMom).setProcess(ProcessType::eDecay);

    children.push_back(std::move(childParticle));
  }
  return children;
}

/** initialize G4RunManager on first call if not done by then */
G4RunManager* 
ActsFatras::Decay::initG4RunManager() const {
  if(G4RunManager::GetRunManager() == nullptr)
  {
	  G4RunManager* runManager = new G4RunManager;
	
	  // initialize here
	  G4VUserPhysicsList *thePL = new QGSP_BERT;

	  runManager->SetUserInitialization(thePL);
	  runManager->SetUserInitialization(new G4DetectorConstruction());

	  // initialize Geant4
	  runManager->Initialize();
	  return runManager;
  } else {
	  return G4RunManager::GetRunManager();
  }
}