// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/EventData/Particle.hpp"
#include "Acts/Utilities/Units.hpp"

#include <TFile.h>
#include <TDirectory.h>
#include <vector>

namespace ActsFatras {

struct Listener {

  const float initialMomentum = 100_GeV;
  
  float invariantMass(const ActsExamples::SimParticle::Vector4& fourVector1,
                    const ActsExamples::SimParticle::Vector4& fourVector2) {
	  ActsExamples::SimParticle::Vector4 sum = fourVector1 + fourVector2;
	  const ActsExamples::SimParticle::Scalar energy = sum[Acts::eEnergy];
	  ActsExamples::SimParticle::Scalar momentum =
		  sum.template segment<3>(Acts::eMom0).norm();
	  return std::sqrt(energy * energy - momentum * momentum);
  }

  std::vector<Particle> operator()(bool soft, const std::vector<Particle>& finalStateParticles, const Particle& interactingParticle, const Particle::Vector4 momentum4) const {
	  // Open the file and prepare everything
	  TFile tf("listener.root", "UPDATE");
	  TDirectory* td;
	  td->cd();
	  auto keys = td->GetListOfKeys();
	  td->cd(keys->First()->GetName());
	  TH1F* histo = nullptr;

	  // Update the nuclear interaction distance
	  td->GetObject("NuclearInteraction", histo);
	  gDirectory->Delete("NuclearInteraction;1");
	  histo->Fill(interactingParticle.pathInL0());
	  gDirectory->WriteObject(histo, "NuclearInteraction");
	  
	  const std::string type = (soft ? "soft" : "hard");
	  td->cd(type.c_str());
	  
	  // Update the multiplicity
	  td->GetObject("Multiplicity", histo);
	  gDirectory->Delete("Multiplicity;1");
	  const unsigned int multiplicity = soft ? finalStateParticles.size() + 1 : finalStateParticles.size();
	  histo->Fill(multiplicity);
	  gDirectory->WriteObject(histo, "Multiplicity");
	  
	  // Update the kinematic plots
	  td->cd(std::to_string(multiplicity).c_str());
	  for(unsigned int i = 0; i < multiplicity; i++)
	  {
		// Index for soft events is shifted by 1
		unsigned int index = soft ? i + 1 : i;
		td->GetObject(("MomentumDistribution_" + std::to_string(index)).c_str(), histo);
		gDirectory->Delete(("MomentumDistribution_" + std::to_string(index) + ";1").c_str());
		histo->Fill(finalStateParticle[i].absMomentum() / initialMomentum;
		gDirectory->WriteObject(histo, ("MomentumDistribution_" + std::to_string(index)).c_str());
		
		td->GetObject(("InvariantMass_" + std::to_string(index)).c_str(), histo);
		gDirectory->Delete(("InvariantMass_" + std::to_string(index) + ";1").c_str());
		histo->Fill(invariantMass(finalStateParticle[i].momentum4(), momentum4));
		gDirectory->WriteObject(histo, ("InvariantMass_" + std::to_string(index)).c_str());
	  }
	  // Set the kinematics of the interacting particle
	  if(soft)
	  {
		td->GetObject("MomentumDistribution_0", histo);
		gDirectory->Delete("MomentumDistribution_0;1");
		histo->Fill(interactingParticle.absMomentum() / initialMomentum;
		gDirectory->WriteObject(histo, "MomentumDistribution_0");
		
		td->GetObject("InvariantMass_0", histo);
		gDirectory->Delete("InvariantMass_0;1");
		histo->Fill(invariantMass(interactingParticle.momentum4(), momentum4));
		gDirectory->WriteObject(histo, "InvariantMass_0");
	  }
	  
	  return finalStateParticles;
  }
  
};
}  // namespace ActsFatras
