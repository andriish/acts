// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include "Acts/Definitions/TrackParametrization.hpp"
#include "TFile.h"
#include "TTree.h"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

struct ParametersAtSurface;
struct TrackSummary;

namespace ActsExamples {

struct Plot {
	
	Plot();
	~Plot();
	
	void mean(const std::vector<TrackSummary>& trackSummaries);
	void scatter(const std::vector<Acts::BoundVector>& localG4Params, 
		const Acts::BoundVector& localPropagatedMean, const std::shared_ptr<const Acts::Surface>& surface);
	
	void plotMean(const std::vector<ParametersAtSurface>& paramAtSurface) const;
	void storeMean(const std::vector<ParametersAtSurface>& paramAtSurface);
	
	Acts::GeometryContext gctx;
	
	TFile* tf{nullptr};
	TTree* tree{nullptr};
	
	std::vector<float> propagatedMean;
	std::vector<float> geant4Mean;
};
}  // namespace ActsExamples
