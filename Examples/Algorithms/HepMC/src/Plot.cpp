// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/HepMC/Plot.hpp"

#include <TGraph.h>
#include <TCanvas.h>
#include "SummaryStruct.hpp"

ActsExamples::Plot::Plot() {
	tf = new TFile("Plot.root", "RECREATE");
	tree = new TTree("Mean", "Mean");
	if (tree == nullptr)
		throw std::bad_alloc();
		
	tree->Branch("propagated", &propagatedMean);
	tree->Branch("geant4", &geant4Mean);
}

ActsExamples::Plot::~Plot() {
	tree->Write();
	tf->Close();
}

void 
ActsExamples::Plot::plotMean(const std::vector<Acts::BoundVector>& props, const std::vector<Acts::BoundVector>& g4s, 
	const std::vector<std::shared_ptr<const Acts::Surface>>& surfaces) const {
	TGraph* tgx = new TGraph();
	TGraph* tgy = new TGraph();
	TGraph* tgphi = new TGraph();
	TGraph* tgtheta = new TGraph();
	TGraph* tgqop = new TGraph();
	TGraph* tgt = new TGraph();
	for(unsigned int i = 0; i < props.size(); i++)
	{
		tgx->SetPoint(i + 1, props[i][Acts::eBoundLoc0], g4s[i][Acts::eBoundLoc0]);
		tgy->SetPoint(i + 1, props[i][Acts::eBoundLoc1], g4s[i][Acts::eBoundLoc1]);
		tgphi->SetPoint(i + 1, props[i][Acts::eBoundPhi], g4s[i][Acts::eBoundPhi]);
		tgtheta->SetPoint(i + 1, props[i][Acts::eBoundTheta], g4s[i][Acts::eBoundTheta]);
		tgqop->SetPoint(i + 1, props[i][Acts::eBoundQOverP], g4s[i][Acts::eBoundQOverP]);
		tgt->SetPoint(i + 1, props[i][Acts::eBoundTime], g4s[i][Acts::eBoundTime]);
	}
	
	TCanvas* tc = new TCanvas();
	tgx->Draw();
	tgy->Draw("same");
	tc->Print("test.png");
}

void 
ActsExamples::Plot::storeMean(const std::vector<Acts::BoundVector>& props, const std::vector<Acts::BoundVector>& g4s,
	const std::vector<std::shared_ptr<const Acts::Surface>>& surfaces) {
	// Store in a TFile
	for(unsigned int i = 0; i < props.size(); i++)
	{
		// Store the surface e.g. for access of radial plots
		const std::shared_ptr<const Acts::Surface> surface = surfaces[i];
		
		// Store the means
		for(unsigned int j = 0; j < Acts::eBoundSize; j++)
		{
			propagatedMean.push_back(props[i][j]);
			geant4Mean.push_back(g4s[i][j]);
		}
		// Write to file
		tree->Fill();
		propagatedMean.clear();
		geant4Mean.clear();
	}
}

void
ActsExamples::Plot::mean(const std::vector<TrackSummary>& trackSummaries) { }
	//~ // Do some plotting
	//~ plotMean(meanProps, meanG4s, surfaces);
	
	//~ // Store the results
	//~ // TODO: move content below if debugged
	//~ storeMean(meanProps, meanG4s, surface);
	

//~ }