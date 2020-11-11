// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <unordered_map>

class G4ParticleDefinition;

namespace ActsFatras
{
  class PDGtoG4Converter
  {
  public:
    /** Default constructor */
    PDGtoG4Converter();
    
/**
   Returns the G4ParticleDefinition of particle with PDG ID pdgCode,
   0 otherwise.
*/
G4ParticleDefinition* getParticleDefinition( int pdgCode) const;

  private:
    /** fills default particles in map of predefined particles */
    void fillPredefinedParticles();
    
/** add a G4 particle to the map of predefined particles */
void addParticle( G4ParticleDefinition* pDef);
    
    /** map from pdg codes to defined Geant4 particles */
    std::unordered_map<int,G4ParticleDefinition*> m_pdgG4ParticleMap;
  };
}