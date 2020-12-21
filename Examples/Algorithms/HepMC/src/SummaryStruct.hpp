// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <vector>
#include "ActsExamples/EventData/SimParticle.hpp"

struct ParametersAtSurface {
	// Mean from propagation, projected onto the surface: Bound and Free representation
	Acts::BoundVector eMeanPropagated;
	Acts::FreeVector eMeanPropagatedFree;
	// Same but for StraightLineStepper
	Acts::BoundVector sMeanPropagated;
	Acts::FreeVector sMeanPropagatedFree;
	
	// Mean from G4 by finding closest points, constructing BoundVectors (by intersecting with the surface and without Cov) and calculating the mean
	Acts::BoundVector meanG4;
	Acts::FreeVector meanG4Free;
	
	// The covariance matrix from the obtained G4 simulations
	Acts::BoundSymMatrix covG4;
	
	// The surface
	std::shared_ptr<const Acts::Surface> surface;
	
	std::vector<Acts::BoundVector> parametersG4;
};

struct TrackSummary {
	// Initial parameters1
	ActsExamples::SimParticle initialParticle;
	
	std::vector<ParametersAtSurface> paramAtSurface;
};