// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Plot.hpp"

#include <TGraph.h>
#include <TCanvas.h>

void 
ActsExamples::Plot::mean(const std::vector<Acts::BoundVector>& meanProp, const std::vector<Acts::BoundVector>& meanG4) {
	TGraph* tgx = new TGraph();
	TGraph* tgy = new TGraph();
	TGraph* tgphi = new TGraph();
	TGraph* tgtheta = new TGraph();
	TGraph* tgqop = new TGraph();
	TGraph* tgt = new TGraph();
	for(unsigned int i = 0; i < meanProp.size(); i++)
	{
		tgx->SetPoint(i + 1, meanProp[i][Acts::eBoundLoc0], meanG4[i][Acts::eBoundLoc0]);
		tgy->SetPoint(i + 1, meanProp[i][Acts::eBoundLoc1], meanG4[i][Acts::eBoundLoc1]);
		tgphi->SetPoint(i + 1, meanProp[i][Acts::eBoundPhi], meanG4[i][Acts::eBoundPhi]);
		tgtheta->SetPoint(i + 1, meanProp[i][Acts::eBoundTheta], meanG4[i][Acts::eBoundTheta]);
		tgqop->SetPoint(i + 1, meanProp[i][Acts::eBoundQOverP], meanG4[i][Acts::eBoundQOverP]);
		tgt->SetPoint(i + 1, meanProp[i][Acts::eBoundTime], meanG4[i][Acts::eBoundTime]);
	}
	
	TCanvas* tc = new TCanvas();
	tgx->Draw();
	tgy->Draw("same");
	tc->Print("test.png");
	
	
	
}