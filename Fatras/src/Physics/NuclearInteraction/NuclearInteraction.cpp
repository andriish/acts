// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Physics/NuclearInteraction/NuclearInteraction.hpp"

namespace ActsFatras {
 
const detail::Parameters& NuclearInteraction::findParameters(double rnd, const detail::Parametrisation& parametrisation, float particleMomentum) const 
{ 
	if(particleMomentum <= parametrisation.begin()->first)
		return parametrisation.begin()->second;
	if(particleMomentum >= parametrisation.end()->first)
		return parametrisation.end()->second;
	
	const auto lowerBound = std::lower_bound(parametrisation.begin(), parametrisation.end(), particleMomentum, 
		[](const std::pair<const float, ActsFatras::detail::Parameters>& params, const float mom){ return params.first < mom; });
	
	const float momentumUpperNeighbour = lowerBound->first;
	const float momentumLowerNeighbour = std::prev(lowerBound, 1)->first;
	const float weight = (momentumUpperNeighbour - particleMomentum) / (momentumUpperNeighbour - momentumLowerNeighbour);
	return (rnd < weight) ? std::prev(lowerBound, 1)->second : lowerBound->second;
}
	
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
	
double NuclearInteraction::sampleContinuousValues(double rnd, const detail::Parameters::CumulativeDistribution& distribution, bool interpolate) const { 	
  // Fast exit
  if(distribution.second.empty()) {
    return 0;
  }
  
  // Find the bin
  const uint32_t int_rnd = UINT32_MAX * rnd;
  const auto it = std::upper_bound(distribution.second.begin(), distribution.second.end(), int_rnd);
  size_t iBin = std::min((size_t) std::distance(distribution.second.begin(), it), distribution.second.size() - 1);
  
  if(interpolate)
  {
	// Interpolate between neighbouring bins and return a diced intermediate value
	const uint32_t basecont = (iBin > 0 ? distribution.second[iBin - 1] : 0);
	const uint32_t dcont = distribution.second[iBin] - basecont;
	return distribution.first[iBin] + (distribution.first[iBin + 1] - distribution.first[iBin]) * (dcont > 0 ? (int_rnd-basecont) / dcont : 0.5);
  }
  else
	return distribution.first[iBin];
}

unsigned int NuclearInteraction::finalStateMultiplicity(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const {
	return sampleDiscreteValues(rnd, distribution);
}

std::pair<ActsFatras::Particle::Scalar, ActsFatras::Particle::Scalar>
NuclearInteraction::globalAngle(ActsFatras::Particle::Scalar phi1, ActsFatras::Particle::Scalar theta1, float phi2, float theta2) const
{	
	const Acts::Vector3F vector2(std::sin(theta2) * std::cos(phi2), std::sin(theta2) * std::sin(phi2), std::cos(theta2));
	Acts::SymMatrix3F rotY = Acts::SymMatrix3F::Zero();
	rotY(0,0) = std::cos(theta1);
	rotY(0,2) = std::sin(theta1);
	rotY(1,1) = 1.;
	rotY(2,0) = -std::sin(theta1);
	rotY(2,2) = std::cos(theta1);
	
	Acts::SymMatrix3F rotZ = Acts::SymMatrix3F::Zero();
	rotZ(0,0) = std::cos(phi1);
	rotZ(0,1) = -std::sin(phi1);
	rotZ(1,0) = std::sin(phi1);
	rotZ(1,1) = std::cos(phi1);
	rotZ(2,2) = 1.;
	
	const Acts::Vector3F vectorSum = rotZ * rotY * vector2;	
	
	const float theta = std::acos(vectorSum.z() / vectorSum.norm());
	const float phi = std::atan2(vectorSum.y(), vectorSum.x());

	return std::make_pair(phi, theta);
}
}  // namespace ActsFatras
