// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Physics/NuclearInteraction/NuclearInteraction.hpp"

namespace ActsFatras {
  
unsigned int NuclearInteraction::sampleDiscreteValues(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const { 
	// Fast exit
  if(distribution.second.empty()) {
    return 0;
  }
  
  // Find the bin
  const uint32_t int_rnd = UINT32_MAX * rnd;
  const auto it = std::upper_bound(distribution.second.begin(), distribution.second.end(), int_rnd);
  size_t iBin = std::min((size_t) std::distance(distribution.second.begin(), it), distribution.second.size() - 1);
  
  // Return the corresponding bin
  return distribution.first[iBin];
  }	
	
double NuclearInteraction::sampleContinuousValues(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const { 	
  // Fast exit
  if(distribution.second.empty()) {
    return 0;
  }
  
  // Find the bin
  const uint32_t int_rnd = UINT32_MAX * rnd;
  const auto it = std::upper_bound(distribution.second.begin(), distribution.second.end(), int_rnd);
  size_t iBin = std::min((size_t) std::distance(distribution.second.begin(), it), distribution.second.size() - 1);
  
  // Interpolate between neighbouring bins and return a diced intermediate value
  const uint32_t basecont = (iBin > 0 ? distribution.second[iBin - 1] : 0);
  const uint32_t dcont = distribution.second[iBin] - basecont;  
  return distribution.first[iBin] + (distribution.first[iBin + 1] - distribution.first[iBin]) * (dcont > 0 ? (int_rnd-basecont) / dcont : 0.5);
  }
}  // namespace ActsFatras
