// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/HepMC/Plot.hpp"

#include <TCanvas.h>
#include "SummaryStruct.hpp"

namespace {
	
Acts::BoundVector
calculateMean(const std::vector<Acts::BoundVector>& positions) {
	Acts::BoundVector mean = Acts::BoundVector::Zero();
	for(const Acts::BoundVector& position : positions)
	{
		mean += position;
	}
	mean /= (double) positions.size();
	return mean;
}
}

ActsExamples::Plot::Plot() {
	tf = TFile::Open("Plot.root", "recreate"); // TODO: open in write mode
	tf->cd();
	
	tgPE1.reserve(10); // nEvents
	tgPE2.reserve(10); // nEvents
	tgPE3.reserve(10); // nEvents
	
	tgPE12->SetLineWidth(0);
	tgPE12->SetMarkerStyle(7);
	tgPE12->SetMarkerColor(kWhite);
	
	tgCum12->SetLineWidth(0);
	tgCum12->SetMarkerStyle(7);
	tgCum12->SetMarkerColor(kWhite);
	
	tgx->SetLineWidth(0);
	tgx->SetMarkerStyle(7);
	tgx->SetMarkerColor(kRed);
	
	tgy->SetLineWidth(0);
	tgy->SetMarkerStyle(7);
	tgy->SetMarkerColor(kWhite);
	
	tgphi->SetLineWidth(0);
	tgphi->SetMarkerStyle(7);
	tgphi->SetMarkerColor(kBlue);
}

ActsExamples::Plot::~Plot() {	
	tf->Close();
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
	tgx->SetMarkerStyle(7);
	tgx->SetMarkerColor(kRed);
	
	tgy->SetLineWidth(0);
	tgy->SetMarkerStyle(8);
	tgy->SetMarkerColor(kBlue);
	
	tgphi->SetLineWidth(0);
	tgphi->SetMarkerStyle(7);
	tgphi->SetMarkerColor(kBlack);
	
	for(unsigned int i = 0; i < paramAtSurface.size(); i++)
	{
		const ParametersAtSurface& pas = paramAtSurface[i];
		const Acts::Vector3D center = pas.surface->center(gctx);
		//~ tgx->SetPoint(i + 1, pas.meanPropagated[Acts::eBoundLoc0], pas.meanG4[Acts::eBoundLoc0] / 100.);
		tgx->SetPoint(i, sqrt(center.x() * center.x() + center.y() * center.y()), pas.eMeanPropagated[Acts::eBoundLoc0] - pas.meanG4[Acts::eBoundLoc0]);
		
		Acts::Vector3D difference = pas.eMeanPropagatedFree.template head<3>(Acts::eFreePos0) - pas.meanG4Free.template head<3>(Acts::eFreePos0);
		tgy->SetPoint(i, sqrt(center.x() * center.x() + center.y() * center.y()), difference.norm() / pas.eMeanPropagatedFree.template head<3>(Acts::eFreePos0).norm());
		
		
		tgphi->SetPoint(i, sqrt(center.x() * center.x() + center.y() * center.y()), difference.norm());
		
		//~ tgy->SetPoint(i + 1, pas.meanPropagated[Acts::eBoundLoc1], pas.meanG4[Acts::eBoundLoc1]);
		//~ tgphi->SetPoint(i, pas.meanPropagated[Acts::eBoundPhi], pas.meanG4[Acts::eBoundPhi]);
		//~ tgtheta->SetPoint(i, pas.meanPropagated[Acts::eBoundTheta], pas.meanG4[Acts::eBoundTheta]);
		//~ tgqop->SetPoint(i, pas.meanPropagated[Acts::eBoundQOverP], pas.meanG4[Acts::eBoundQOverP]);
		//~ tgt->SetPoint(i, pas.meanPropagated[Acts::eBoundTime], pas.meanG4[Acts::eBoundTime]);
	}
	
	tgx->SetTitle(";Transversal distance of surface center [mm];#Delta Loc0 [mm]");
	
	TCanvas* tc = new TCanvas();
	tgx->Draw();
	tgy->Draw("psame");
	tgphi->Draw("psame");
	tc->Print("test.png");
}

void
ActsExamples::Plot::plotCumulativeMean(const std::vector<TrackSummary>& trackSummaries) {
	
	TGraph* tg1 = new TGraph();
	TGraph* tg2 = new TGraph();
	
	tg1->SetLineColor(kRed - tgPE1.size());
	tg1->SetMarkerStyle(8);
	tg1->SetMarkerColor(kRed - tgPE1.size());
	
	tg2->SetLineColor(kBlue - tgPE2.size());
	tg2->SetMarkerStyle(8);
	tg2->SetMarkerColor(kBlue - tgPE2.size());
	
	for(const auto& ts : trackSummaries)
	{
		const auto& paramAtSurface = ts.paramAtSurface;
		for(unsigned int i = 0; i < paramAtSurface.size(); i++)
		{		
			const ParametersAtSurface& pas = paramAtSurface[i];
			const Acts::Vector3D center = pas.surface->center(gctx);
			const double radius = sqrt(center.x() * center.x() + center.y() * center.y());
		
			parametersG4Cumulative[pas.surface].insert(parametersG4Cumulative[pas.surface].end(), pas.parametersG4.begin(), pas.parametersG4.end()); 
			const Acts::BoundVector meanG4 = calculateMean(parametersG4Cumulative[pas.surface]);
			
			tgCum12->SetPoint(2. * (i + nPoints), radius, pas.eMeanPropagated[Acts::eBoundLoc0] - meanG4[Acts::eBoundLoc0]);
			tg1->SetPoint(i, radius, pas.eMeanPropagated[Acts::eBoundLoc0] - meanG4[Acts::eBoundLoc0]);
			
			const Acts::FreeVector meanG4Free = Acts::detail::transformBoundToFreeParameters(*pas.surface, gctx, meanG4);
			Acts::Vector3D difference = pas.eMeanPropagatedFree.template head<3>(Acts::eFreePos0) - meanG4Free.template head<3>(Acts::eFreePos0);
			
			tgCum12->SetPoint(2. * (i + nPoints) + 1, radius, difference.norm());
			tg2->SetPoint(i, radius, difference.norm());
		}
		 nPoints += paramAtSurface.size();	
	}
	
	tgCum1.push_back(tg1);
	tgCum2.push_back(tg2);
	
	TCanvas* tc = new TCanvas();
	tgCum12->SetTitle(";Transversal distance of surface center [mm];#Delta Loc0 [mm]");
	tgCum12->Draw();
	for(TGraph* t : tgCum1)
		t->Draw("plsame");
	tc->Print("testCum.png");
	delete(tc);
	
	tc = new TCanvas();
	tgCum12->SetTitle(";Transversal distance of surface center [mm];#Delta global position [mm]");
	tgCum12->Draw();
	for(TGraph* t : tgCum2)
		t->Draw("plsame");
	tc->Print("testCum2.png");
	delete(tc);	
}

void
ActsExamples::Plot::mean(const std::vector<TrackSummary>& trackSummaries) {
	
	TGraph* tg1 = new TGraph();
	//~ TGraph* tg1 = tgPE1.empty() ? new TGraph() : new TGraph(*tgPE1.back());
	//~ if(!tgPE1.empty())
	//~ {
		//~ for(int i = 0; i < tgPE1.back()->GetN());
			//~ tg1->SetPoint(i, tgPE1.back()->GetPoint
	//~ }
	TGraph* tg2 = new TGraph();
	TGraph* tg3 = new TGraph();
	
	tg1->SetLineColor(kRed - tgPE1.size());
	tg1->SetMarkerStyle(8);
	tg1->SetMarkerColor(kRed - tgPE1.size());
	
	tg2->SetLineColor(kBlue - tgPE2.size());
	tg2->SetMarkerStyle(8);
	tg2->SetMarkerColor(kBlue - tgPE2.size());
	
	tg3->SetLineColor(kGreen - tgPE3.size());
	tg3->SetMarkerStyle(8);
	tg3->SetMarkerColor(kGreen - tgPE3.size());
	
	for(const auto& ts : trackSummaries)
	{		
		const auto& paramAtSurface = ts.paramAtSurface;
		for(unsigned int i = 0; i < paramAtSurface.size(); i++)
		{
			const ParametersAtSurface& pas = paramAtSurface[i];
			const Acts::Vector3D center = pas.surface->center(gctx);
			Acts::Vector3D difference = pas.eMeanPropagatedFree.template head<3>(Acts::eFreePos0) - pas.meanG4Free.template head<3>(Acts::eFreePos0);
			const double radius = sqrt(center.x() * center.x() + center.y() * center.y());
			
			tgPE12->SetPoint(2. * (i + nPoints), radius, pas.eMeanPropagated[Acts::eBoundLoc0] - pas.meanG4[Acts::eBoundLoc0]);
			tg1->SetPoint(i, radius, pas.eMeanPropagated[Acts::eBoundLoc0] - pas.meanG4[Acts::eBoundLoc0]);
			
			tg3->SetPoint(i, radius, pas.sMeanPropagated[Acts::eBoundLoc0] - pas.meanG4[Acts::eBoundLoc0]);
			
			tgPE12->SetPoint(2. * (i + nPoints) + 1, radius, difference.norm());
			tg2->SetPoint(i, radius, difference.norm());
		}
		 nPoints += paramAtSurface.size();
	}
	
	tgPE1.push_back(tg1);
	tgPE2.push_back(tg2);
	tgPE3.push_back(tg3);
	
	TCanvas* tc = new TCanvas();
	tgPE12->SetTitle(";Transversal distance of surface center [mm];#Delta Loc0 [mm]");
	tgPE12->Draw();
	for(TGraph* t : tgPE1)
		t->Draw("plsame");
	tc->Print("test.png");
	delete(tc);
	
	tc = new TCanvas();
	tgPE12->Draw();
	for(TGraph* t : tgPE3)
		t->Draw("plsame");
	tc->Print("test1.png");
	delete(tc);
	
	tc = new TCanvas();
	tgPE12->SetTitle(";Transversal distance of surface center [mm];#Delta global position [mm]");
	tgPE12->Draw();
	for(TGraph* t : tgPE2)
		t->Draw("plsame");
	tc->Print("test2.png");
	delete(tc);
	
	plotCumulativeMean(trackSummaries);
}

void 
ActsExamples::Plot::scatter(const std::vector<Acts::BoundVector>& localG4Params, 
		const Acts::BoundVector& eMeanPropagated, const Acts::BoundVector& sMeanPropagated, const std::shared_ptr<const Acts::Surface>& surface) {
	
	Acts::BoundVector mean2 = Acts::BoundVector::Zero();
	for(const Acts::BoundVector& position : localG4Params)
	{
		mean2 += position;
	}
	mean2 /= (double) localG4Params.size();
	
	std::vector<std::vector<float>> parametersHit; // [parameter][event]
	parametersHit.resize(Acts::eBoundSize);
	std::vector<std::vector<float>> parametersMiss; // [parameter][event]
	parametersMiss.resize(Acts::eBoundSize);
	std::vector<float> eMean;
	eMean.resize(Acts::eBoundSize);
	std::vector<float> sMean;
	sMean.resize(Acts::eBoundSize);
	
	for(unsigned int i = 0; i < Acts::eBoundSize; i++)
	{
		eMean[i] = eMeanPropagated[i];
		sMean[i] = sMeanPropagated[i];
	}
	
	const Acts::Vector3D center = surface->center(gctx);

	tf->WriteObject(&eMean, ("MeanES_" + std::to_string((int) center.x()) + "_" 
		+ std::to_string((int) center.y()) + "_" + std::to_string((int) center.z())).c_str());
	tf->WriteObject(&sMean, ("MeanSLS_" + std::to_string((int) center.x()) + "_" 
		+ std::to_string((int) center.y()) + "_" + std::to_string((int) center.z())).c_str());
	
  
	TGraph* hit = new TGraph();
	TGraph* miss = new TGraph();
	TGraph* both = new TGraph();
	TGraph* meanBoth = new TGraph(1);
	TGraph* eMeanProp = new TGraph(1);
	TGraph* sMeanProp = new TGraph(1);
	
	hit->SetLineWidth(0);
	hit->SetMarkerStyle(7);
	hit->SetMarkerColor(kBlue);
	
	miss->SetLineWidth(0);
	miss->SetMarkerStyle(7);
	miss->SetMarkerColor(kRed);
	
	both->SetLineWidth(0);
	both->SetMarkerStyle(7);
	both->SetMarkerColor(kWhite);
	//~ both->SetTitle(";Loc0 [mm];Loc1 [mm]");
	
	meanBoth->SetLineWidth(0);
	meanBoth->SetMarkerStyle(8);
	meanBoth->SetMarkerColor(kGreen);
	
	eMeanProp->SetLineWidth(0);
	eMeanProp->SetMarkerStyle(8);
	eMeanProp->SetMarkerColor(kBlack);
	
	sMeanProp->SetLineWidth(0);
	sMeanProp->SetMarkerStyle(8);
	sMeanProp->SetMarkerColor(kBlack);
	
	unsigned int h = 0;
	unsigned int m = 0;
	
	for(const auto& p : localG4Params)
	{
		both->SetPoint(h + m, p[Acts::eBoundLoc0], p[Acts::eBoundLoc1]);
		if(surface->insideBounds(p.template head<2>()))
		{
			hit->SetPoint(h, p[Acts::eBoundLoc0], p[Acts::eBoundLoc1]);
			h++;
			for(unsigned int i = 0; i < Acts::eBoundSize; i++)
				parametersHit[i].push_back(p[i]);
		}
		else
		{
			miss->SetPoint(m, p[Acts::eBoundLoc0], p[Acts::eBoundLoc1]);
			m++;
			for(unsigned int i = 0; i < Acts::eBoundSize; i++)
				parametersMiss[i].push_back(p[i]);
		}
	}
	
	meanBoth->SetPoint(0, mean2[Acts::eBoundLoc0], mean2[Acts::eBoundLoc1]);
	eMeanProp->SetPoint(0, eMeanPropagated[Acts::eBoundLoc0], eMeanPropagated[Acts::eBoundLoc1]);
	sMeanProp->SetPoint(0, sMeanPropagated[Acts::eBoundLoc0], sMeanPropagated[Acts::eBoundLoc1]);
	

	TCanvas* tc = new TCanvas();
	both->Draw();
	if(m > 0)
		miss->Draw("psame");
	if(h > 0)
		hit->Draw("psame");
	meanBoth->Draw("psame");
	eMeanProp->Draw("psame");
	tc->Print(("testScatter" + std::to_string((int) center.x()) + "_" 
		+ std::to_string((int) center.y()) + "_" + std::to_string((int) center.z()) + ".png").c_str());
	
	std::vector<double> centerVector{center.x(), center.y(), center.z()};
	tf->WriteObject(&centerVector, ("SurfaceCenter_" + std::to_string((int) center.x()) + "_" 
		+ std::to_string((int) center.y()) + "_" + std::to_string((int) center.z())).c_str());
	
	for(unsigned int i = 0; i < Acts::eBoundSize; i++)
	{
		tf->WriteObject(&parametersHit[i], ("Hit_" + std::to_string(i) + "_" + std::to_string((int) center.x()) + "_" 
			+ std::to_string((int) center.y()) + "_" + std::to_string((int) center.z())).c_str());
		tf->WriteObject(&parametersMiss[i], ("Miss_" + std::to_string(i) + "_" + std::to_string((int) center.x()) + "_" 
			+ std::to_string((int) center.y()) + "_" + std::to_string((int) center.z())).c_str());
	}
}