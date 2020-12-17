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
	tf = new TFile("Plot.root", "recreate"); // TODO: open in write mode
	tf->cd();
	tree = new TTree("Mean", "Mean");
	if (tree == nullptr)
		throw std::bad_alloc();
		
	tree->Branch("propagated", &propagatedMean);
	tree->Branch("geant4", &geant4Mean);
}

ActsExamples::Plot::~Plot() {
	//~ tf->cd();
	//~ tree->Write();
	//~ tf->Close();
}

void 
ActsExamples::Plot::plotMean(const std::vector<ParametersAtSurface>& paramAtSurface) const {
	TGraph* tgx = new TGraph();
	TGraph* tgy = new TGraph();
	TGraph* tgphi = new TGraph();
	TGraph* tgtheta = new TGraph();
	TGraph* tgqop = new TGraph();
	TGraph* tgt = new TGraph();

	tgx->SetLineWidth(0);
	tgx->SetMarkerStyle(8);
	tgx->SetMarkerColor(kRed);
	
	tgy->SetLineWidth(0);
	tgy->SetMarkerStyle(8);
	tgy->SetMarkerColor(kBlue);
	
	tgphi->SetLineWidth(0);
	tgphi->SetMarkerStyle(8);
	tgphi->SetMarkerColor(kBlack);
	
	for(unsigned int i = 0; i < paramAtSurface.size(); i++)
	{
		const ParametersAtSurface& pas = paramAtSurface[i];
		const Acts::Vector3D center = pas.surface->center(gctx);
		//~ tgx->SetPoint(i + 1, pas.meanPropagated[Acts::eBoundLoc0], pas.meanG4[Acts::eBoundLoc0] / 100.);
		tgx->SetPoint(i, sqrt(center.x() * center.x() + center.y() * center.y()), pas.meanPropagated[Acts::eBoundLoc0] - pas.meanG4[Acts::eBoundLoc0]);
		
		Acts::Vector3D difference = pas.meanPropagatedFree.template head<3>(Acts::eFreePos0) - pas.meanG4Free.template head<3>(Acts::eFreePos0);
		tgy->SetPoint(i, sqrt(center.x() * center.x() + center.y() * center.y()), difference.norm() / pas.meanPropagatedFree.template head<3>(Acts::eFreePos0).norm());
		
		
		tgphi->SetPoint(i, sqrt(center.x() * center.x() + center.y() * center.y()), difference.norm());
		
		//~ tgy->SetPoint(i + 1, pas.meanPropagated[Acts::eBoundLoc1], pas.meanG4[Acts::eBoundLoc1]);
		//~ tgphi->SetPoint(i, pas.meanPropagated[Acts::eBoundPhi], pas.meanG4[Acts::eBoundPhi]);
		tgtheta->SetPoint(i, pas.meanPropagated[Acts::eBoundTheta], pas.meanG4[Acts::eBoundTheta]);
		tgqop->SetPoint(i, pas.meanPropagated[Acts::eBoundQOverP], pas.meanG4[Acts::eBoundQOverP]);
		tgt->SetPoint(i, pas.meanPropagated[Acts::eBoundTime], pas.meanG4[Acts::eBoundTime]);
	}
	
	tgx->SetTitle(";Transversal distance of surface center [mm];#Delta Loc0 [mm]");
	
	TCanvas* tc = new TCanvas();
	tgx->Draw();
	tgy->Draw("psame");
	tgphi->Draw("psame");
	tc->Print("test.png");
}

void 
ActsExamples::Plot::storeMean(const std::vector<ParametersAtSurface>& paramAtSurface) {
	// Store in a TFile
	for(const ParametersAtSurface& pas : paramAtSurface)
	{
		// Store the surface e.g. for access of radial plots
		//~ const std::shared_ptr<const Acts::Surface> surface = surfaces[i];
		
		// Store the means
		for(unsigned int j = 0; j < Acts::eBoundSize; j++)
		{
			propagatedMean.push_back(pas.meanPropagated[j]);
			geant4Mean.push_back(pas.meanG4[j]);
		}
		// Write to file
		tree->Fill();
		propagatedMean.clear();
		geant4Mean.clear();
	}
}

void
ActsExamples::Plot::mean(const std::vector<TrackSummary>& trackSummaries) {
	
	TGraph* tgx = new TGraph();
	TGraph* tgy = new TGraph();
	TGraph* tgphi = new TGraph();
	TGraph* tgtheta = new TGraph();
	TGraph* tgqop = new TGraph();
	TGraph* tgt = new TGraph();

	tgx->SetLineWidth(0);
	tgx->SetMarkerStyle(8);
	tgx->SetMarkerColor(kRed);
	
	tgy->SetLineWidth(0);
	tgy->SetMarkerStyle(8);
	tgy->SetMarkerColor(kBlue);
	
	tgphi->SetLineWidth(0);
	tgphi->SetMarkerStyle(8);
	tgphi->SetMarkerColor(kBlack);
	
	for(const auto& ts : trackSummaries)
	{		
		const auto& paramAtSurface = ts.paramAtSurface;
		for(unsigned int i = 0; i < paramAtSurface.size(); i++)
		{
			const ParametersAtSurface& pas = paramAtSurface[i];
			const Acts::Vector3D center = pas.surface->center(gctx);
			//~ tgx->SetPoint(i + 1, pas.meanPropagated[Acts::eBoundLoc0], pas.meanG4[Acts::eBoundLoc0] / 100.);
			tgx->SetPoint(i, sqrt(center.x() * center.x() + center.y() * center.y()), pas.meanPropagated[Acts::eBoundLoc0] - pas.meanG4[Acts::eBoundLoc0]);
			
			Acts::Vector3D difference = pas.meanPropagatedFree.template head<3>(Acts::eFreePos0) - pas.meanG4Free.template head<3>(Acts::eFreePos0);
			tgy->SetPoint(i, sqrt(center.x() * center.x() + center.y() * center.y()), difference.norm() / pas.meanPropagatedFree.template head<3>(Acts::eFreePos0).norm());
			
			
			tgphi->SetPoint(i, sqrt(center.x() * center.x() + center.y() * center.y()), difference.norm());
			
			//~ tgy->SetPoint(i + 1, pas.meanPropagated[Acts::eBoundLoc1], pas.meanG4[Acts::eBoundLoc1]);
			//~ tgphi->SetPoint(i, pas.meanPropagated[Acts::eBoundPhi], pas.meanG4[Acts::eBoundPhi]);
			tgtheta->SetPoint(i, pas.meanPropagated[Acts::eBoundTheta], pas.meanG4[Acts::eBoundTheta]);
			tgqop->SetPoint(i, pas.meanPropagated[Acts::eBoundQOverP], pas.meanG4[Acts::eBoundQOverP]);
			tgt->SetPoint(i, pas.meanPropagated[Acts::eBoundTime], pas.meanG4[Acts::eBoundTime]);
		}	
	}
	tgx->SetTitle(";Transversal distance of surface center [mm];#Delta Loc0 [mm]");
	
	TCanvas* tc = new TCanvas();
	tgx->Draw();
	tgy->Draw("psame");
	tgphi->Draw("psame");
	tc->Print("test.png");
}

void 
ActsExamples::Plot::scatter(const std::vector<Acts::BoundVector>& localG4Params, 
		const Acts::BoundVector& localPropagatedMean, const std::shared_ptr<const Acts::Surface>& surface) {
	TGraph* hit = new TGraph();
	TGraph* miss = new TGraph();
	TGraph* both = new TGraph();
	TGraph* meanBoth = new TGraph(1);
	TGraph* meanHit = new TGraph(1);
	TGraph* meanProp = new TGraph(1);
	
	hit->SetLineWidth(0);
	hit->SetMarkerStyle(7);
	hit->SetMarkerColor(kBlue);
	
	miss->SetLineWidth(0);
	miss->SetMarkerStyle(7);
	miss->SetMarkerColor(kRed);
	
	both->SetLineWidth(0);
	both->SetMarkerStyle(7);
	both->SetMarkerColor(kWhite);
	both->SetTitle(";Loc0 [mm];Loc1 [mm]");
	
	meanBoth->SetLineWidth(0);
	meanBoth->SetMarkerStyle(8);
	meanBoth->SetMarkerColor(kGreen);
	
	meanHit->SetLineWidth(0);
	meanHit->SetMarkerStyle(8);
	meanHit->SetMarkerColor(kMagenta);
	
	meanProp->SetLineWidth(0);
	meanProp->SetMarkerStyle(8);
	meanProp->SetMarkerColor(kBlack);
	
	unsigned int h = 0;
	unsigned int m = 0;
	
	Acts::BoundVector mean2 = Acts::BoundVector::Zero();
	for(const Acts::BoundVector& position : localG4Params)
	{
		mean2 += position;
	}
	mean2 /= (double) localG4Params.size();
	
	Acts::BoundVector mean1 = Acts::BoundVector::Zero();
	for(const auto& p : localG4Params)
	{
		both->SetPoint(h + m, p[Acts::eBoundLoc0], p[Acts::eBoundLoc1]);
		if(surface->insideBounds(p.template head<2>()))
		{
			hit->SetPoint(h, p[Acts::eBoundLoc0], p[Acts::eBoundLoc1]);
			h++;
			mean1 += p;
		}
		else
		{
			miss->SetPoint(m, p[Acts::eBoundLoc0], p[Acts::eBoundLoc1]);
			m++;
		}
	}
	mean1 /= (double) (h + 1);
	
	meanBoth->SetPoint(0, mean2[Acts::eBoundLoc0], mean2[Acts::eBoundLoc1]);
	meanHit->SetPoint(0, mean1[Acts::eBoundLoc0], mean1[Acts::eBoundLoc1]);
	meanProp->SetPoint(0, localPropagatedMean[Acts::eBoundLoc0], localPropagatedMean[Acts::eBoundLoc1]);
	
	const Acts::Vector3D center = surface->center(gctx);
	//~ const double radius = sqrt(center.x() * center.x() + center.y() * center.y());
std::cout << "Center: " << center.transpose() << std::endl;
std::cout << "Mean: " << mean2[Acts::eBoundLoc0] << ", " << mean2[Acts::eBoundLoc1] 
<< " | " << localPropagatedMean[Acts::eBoundLoc0] << ", " << localPropagatedMean[Acts::eBoundLoc1] << std::endl;
	TCanvas* tc = new TCanvas();
	both->Draw();
	if(m > 0)
		miss->Draw("psame");
	if(h > 0)
		hit->Draw("psame");
	meanBoth->Draw("psame");
	meanHit->Draw("psame");
	meanProp->Draw("psame");
	tc->Print(("testScatter" + std::to_string(center.x()) + "_" + std::to_string(center.y()) + "_" + std::to_string(center.z()) + ".png").c_str());
	
	tf->Write();
}