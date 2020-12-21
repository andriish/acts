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
#include <TGraph.h>
#include <map>

struct ParametersAtSurface;
struct TrackSummary;

namespace ActsExamples {

struct Plot {
	
	Plot();
	~Plot();
	
	void mean(const std::vector<TrackSummary>& trackSummaries);
	void scatter(const std::vector<Acts::BoundVector>& localG4Params, 
		const Acts::BoundVector& eMeanPropagated, const Acts::BoundVector& sMeanPropagated, const std::shared_ptr<const Acts::Surface>& surface);
	
	void plotMean(const std::vector<ParametersAtSurface>& paramAtSurface) const;
	void plotCumulativeMean(const std::vector<TrackSummary>& trackSummaries);
	
	Acts::GeometryContext gctx;
	
	TFile* tf{nullptr};
	
	std::vector<TGraph*> tgPE1; // PE = per event
	std::vector<TGraph*> tgPE2;
	std::vector<TGraph*> tgPE3;
	TGraph* tgPE12 = new TGraph();
	
	std::map<std::shared_ptr<const Acts::Surface>, std::vector<Acts::BoundVector>> parametersG4Cumulative;
	std::vector<TGraph*> tgCum1;
	std::vector<TGraph*> tgCum2;
	TGraph* tgCum12 = new TGraph();
	
	TGraph* tgx = new TGraph(); // tgxDiff
	TGraph* tgy = new TGraph(); // tgDiffGlobalPosition
	TGraph* tgphi = new TGraph(); // tgRelativeDiffGlobalPosition
	TGraph* tgtheta = new TGraph();
	TGraph* tgqop = new TGraph();
	TGraph* tgt = new TGraph();
	unsigned int nPoints = 0;

};
}  // namespace ActsExamples
