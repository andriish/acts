// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

//~ // CLHEP
//~ #include "CLHEP/Units/SystemOfUnits.h"
//~ #include "CLHEP/Units/PhysicalConstants.h"

#include <vector>
#include <cmath>
#include <random>
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Physics/Decay/PDGtoG4Converter.hpp"
#include "G4ParticleDefinition.hh"

class G4RunManager;

namespace ActsFatras {

  class Decay {
    public:
      
      Decay();
      
/** free path estimator (-1 for stable particle) */
template <typename generator_t>
double freeLifeTime(generator_t& generator, const Particle& isp) const;

/** decay handling secondaries */
template <typename generator_t>
std::vector<Particle> operator()(generator_t& generator, const Acts::MaterialSlab& /*slab*/, const Particle& isp) const;

/** decay */
std::vector<Particle> decayParticle(const Particle& parent) const;

   private:
/** initialize G4RunManager on first call if not done by then */
G4RunManager* initG4RunManager() const;

      mutable G4RunManager*                m_g4RunManager;         //!< for dummy G4 initialization               
                      
      PDGtoG4Converter         m_pdgToG4Conv;           //!< Handle for the  PDGToG4Particle converter tool

   };               
                    
}

template <typename generator_t>
double 
ActsFatras::Decay::freeLifeTime(generator_t& generator, const ActsFatras::Particle& isp) const {
  // get the particle properties    
  const int pdgCode =  isp.pdg();
  
  if(std::abs(pdgCode) == 13) return -1.;
  
  G4ParticleDefinition* pDef = m_pdgToG4Conv.getParticleDefinition(pdgCode);
  
  if(!pDef || pDef->GetPDGStable()) {
    return -1.;
  }
    
  // take momentum from ParticleState rather than associated truth
  Particle::Vector4 particleMom = isp.momentum4();

  // get average lifetime
  constexpr double convertTime = Acts::UnitConstants::mm / CLHEP::s;
  const double tau = pDef->GetPDGLifeTime() * convertTime;
  // sample the lifetime
  std::uniform_real_distribution<double> uniformDistribution {0., 1.};
  const double lifeTime = -tau * log(uniformDistribution(generator));
  
  return lifeTime;
}

template <typename generator_t>
std::vector<ActsFatras::Particle> 
ActsFatras::Decay::operator()(generator_t& generator, const Acts::MaterialSlab& /*slab*/, const ActsFatras::Particle& isp) const{

  // perform the decay 
  const std::vector<Particle> decayProducts = decayParticle(isp);
  
	// TODO: free path must be set for decay products
	// TODO: initial particle must die
	
	return decayProducts;
}