// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

//~ class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;

// Geant4
#include "G4VUserDetectorConstruction.hh"

namespace ActsFatras {
class G4DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  ~G4DetectorConstruction();

  G4VPhysicalVolume* Construct();

private:
  void dummyDetector();

  // Logical volume
  G4LogicalVolume*          m_worldLog = nullptr;

  // Physical volume
  G4VPhysicalVolume*        m_worldPhys = nullptr;
};
}