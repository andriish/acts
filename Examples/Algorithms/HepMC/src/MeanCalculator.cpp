// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"

#include "ActsExamples/HepMC/MeanCalculator.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Particle.hpp"

#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

#include "Acts/Propagator/Navigator.hpp"

#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include "SummaryStruct.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"

namespace {

std::vector<ActsExamples::SimParticle>
collectG4Steps(const HepMC3::GenEvent& event, int trackID) {
	std::vector<ActsExamples::SimParticle> g4Steps;
	for (const auto& vertex : event.vertices()) {
		const auto attribute = vertex->attribute<HepMC3::StringAttribute>("Material");
		if(attribute == nullptr)
		{
			continue;
		}
		const auto material = attribute->value();
		if(material == "NoMaterial" || material == "Vacuum" || material == "Air")
			continue;
		for(const auto& particle :vertex->particles_out()) {
			if(particle->attribute<HepMC3::IntAttribute>("TrackID")->value() == trackID)
			{
				const auto& posVtx = vertex->position();
				const Acts::Vector4D position(posVtx.x(), posVtx.y(), posVtx.z(), posVtx.t());
				
				const auto& momPart = particle->momentum();
				const Acts::Vector3D momentum(momPart.x(), momPart.y(), momPart.z());
				
				ActsExamples::SimParticle g4Particle;
				g4Particle.setPosition4(position).setDirection(momentum.normalized()).setAbsoluteMomentum(momentum.norm());
				g4Steps.push_back(g4Particle);
				break;
			}
		}
	}
	return g4Steps;
}

std::optional<Acts::BoundVector>
findClosestPoint(const std::vector<ActsExamples::SimParticle>& g4Steps, const std::shared_ptr<const Acts::Surface> surface, const Acts::GeometryContext& gctx) {
	std::vector<std::pair<double, Acts::BoundVector>> pathLengthPosition;
	Acts::BoundaryCheck bCheck(false);
	for(const ActsExamples::SimParticle& g4Step : g4Steps)
	{
		const Acts::SurfaceIntersection intersection = surface->intersect(gctx, g4Step.position(), g4Step.unitDirection(), bCheck);
		if(intersection)
		{
			const Acts::Vector3D pos3 = intersection.intersection.position;
			Acts::FreeVector freeParams;
			freeParams[Acts::eFreePos0] = pos3[Acts::eX];
			freeParams[Acts::eFreePos1] = pos3[Acts::eY];
			freeParams[Acts::eFreePos2] = pos3[Acts::eZ];
			freeParams[Acts::eFreeTime] = g4Step.time();
			freeParams[Acts::eFreeDir0] = g4Step.unitDirection()[Acts::eMom0];
			freeParams[Acts::eFreeDir1] = g4Step.unitDirection()[Acts::eMom1];
			freeParams[Acts::eFreeDir2] = g4Step.unitDirection()[Acts::eMom2];
			freeParams[Acts::eFreeQOverP] = (g4Step.charge() == 0. ? 1. : g4Step.charge()) / g4Step.absoluteMomentum();
			
			Acts::BoundVector params = Acts::detail::transformFreeToBoundParameters(freeParams, *surface, gctx);
		
			const double pathLength = intersection.intersection.pathLength;
			pathLengthPosition.push_back(std::make_pair(pathLength, params));
		}
	}
	if(!pathLengthPosition.empty())
	{
		const auto closest = std::min_element(pathLengthPosition.begin(), pathLengthPosition.end(), 
				[&](const std::pair<double, Acts::BoundVector>& pos1, const std::pair<double, Acts::BoundVector>& pos2) 
					{ return  std::abs(pos1.first) < std::abs(pos2.first); });

		return closest->second;
	} else
		return std::nullopt;
}

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

Acts::BoundSymMatrix
calculateCovariance(const std::vector<Acts::BoundVector>& params, const Acts::BoundVector& mean) {
  // Calculate the covariance matrix
  Acts::BoundSymMatrix covariance = Acts::BoundSymMatrix::Zero(Acts::eBoundSize, Acts::eBoundSize);
  for (unsigned int i = 0; i < Acts::eBoundSize; i++)
    for (unsigned int j = 0; j < Acts::eBoundSize; j++)
      for (unsigned int k = 0; k < params.size(); k++)
        covariance(i, j) += (params[k][i] - mean[i]) * (params[k][j] - mean[j]);
  covariance /= (float)params.size();
  
  return covariance;
}

std::pair<Acts::BoundVector, Acts::FreeVector>
convertStepToParameters(Acts::GeometryContext& gctx, const Acts::detail::Step& step, double charge) {
	// Calculate the value of the mean on the surface
	Acts::Vector3D dir = step.momentum.normalized();
	Acts::FreeVector freeProp;
	freeProp[Acts::eFreePos0] = step.position[Acts::eX];
	freeProp[Acts::eFreePos1] = step.position[Acts::eY];
	freeProp[Acts::eFreePos2] = step.position[Acts::eZ];
	freeProp[Acts::eFreeTime] = step.time;
	freeProp[Acts::eFreeDir0] = dir[Acts::eMom0];
	freeProp[Acts::eFreeDir1] = dir[Acts::eMom1];
	freeProp[Acts::eFreeDir2] = dir[Acts::eMom2];
	freeProp[Acts::eFreeQOverP] = (charge == 0. ? 1. : charge) / step.momentum.norm();
	// Build the local parameters
	const Acts::BoundVector localPropagatedMean = Acts::detail::transformFreeToBoundParameters(freeProp, *step.surface, gctx);
	
	return std::make_pair(localPropagatedMean, freeProp);
}

}  // namespace

ActsExamples::MeanCalculator::~MeanCalculator() {}

ActsExamples::MeanCalculator::MeanCalculator(
    ActsExamples::MeanCalculator::Config&& cnf,
    Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("MeanCalculator", level),
      m_cfg(std::move(cnf)) {
  if (m_cfg.inputEvents.empty()) {
    throw std::invalid_argument("Missing input event collection");
  }
  if (m_cfg.inputParticles.empty()) {
	throw std::invalid_argument("Missing input event collection");
  }
}

ActsExamples::ProcessCode ActsExamples::MeanCalculator::execute(
    const ActsExamples::AlgorithmContext& context) const {
  // Retrieve the initial particles
  const auto events =
      context.eventStore.get<std::vector<HepMC3::GenEvent>>(m_cfg.inputEvents);

  // Retrieve the initial particles
  const auto initialParticles =
      context.eventStore.get<ActsExamples::SimParticleContainer>(
          m_cfg.inputParticles);
  
  std::vector<TrackSummary> summaries;
  
  // The stepper
  Acts::NullBField bfield;  
  Acts::EigenStepper eStepper(bfield);
  Acts::StraightLineStepper slStepper;
  
  // The Navigator
  Acts::Navigator navigator(m_cfg.trackingGeometry);
  
  // The Propagator
  Acts::Propagator ePropagator(eStepper, navigator);
  Acts::Propagator slPropagator(slStepper, navigator);
  
  // The options
  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::PropagatorOptions<Acts::ActionList<Acts::detail::SteppingLogger>, Acts::AbortList<Acts::EndOfWorldReached>> options(gctx, mctx, Acts::getDummyLogger());
  options.maxStepSize = m_cfg.maxStepSize;
  options.maxSteps = m_cfg.maxSteps;
  options.tolerance = m_cfg.tolerance;
  
  // Loop over initial particles
  for(const ActsExamples::SimParticle& initialParticle : initialParticles)
  {
	    // Storage of the summary of a single track
	    TrackSummary summary;
	    summary.initialParticle = initialParticle;
	  
	    // Propagate the mean
		Acts::CurvilinearTrackParameters mean(initialParticle.fourPosition(), 
			initialParticle.unitDirection(), initialParticle.charge(), initialParticle.absoluteMomentum());
		const auto& eResult = ePropagator.propagate(mean, options).value(); //result.ok()
		const auto eStepperLog = eResult.get<typename Acts::detail::SteppingLogger::result_type>();
		const auto& slResult = slPropagator.propagate(mean, options).value();
		const auto slStepperLog = slResult.get<typename Acts::detail::SteppingLogger::result_type>();
				
		/// Find the surfaces first
		// Walk over each step
		for(const auto& step : eStepperLog.steps)
		{
			// Only care about surfaces
			if(!step.surface)
				continue;
			const Acts::Vector3D surfaceCenter = step.surface->center(gctx);
			if(sqrt(surfaceCenter.x() * surfaceCenter.x() + surfaceCenter.y() * surfaceCenter.y())  == 0.)
				continue;

			ParametersAtSurface surfaceParameters;
			
			// Get the corresponding mean from the straight line stepper
			for(const auto& slStep : slStepperLog.steps)
			{
				if(step.surface == slStep.surface)
				{
					const auto params = convertStepToParameters(gctx, slStep, mean.charge());
					surfaceParameters.sMeanPropagated = std::move(params.first);
					surfaceParameters.sMeanPropagatedFree = std::move(params.second);
				}
			}


			// Store the corresponding for surface
			surfaceParameters.surface = step.surface;

			// Store the EigenStepper parameters at the surface
			const auto params = convertStepToParameters(gctx, step, mean.charge());
			surfaceParameters.eMeanPropagated = std::move(params.first);
			surfaceParameters.eMeanPropagatedFree = std::move(params.second);
			
			// Now find the corresponding G4 steps
			std::vector<Acts::BoundVector> localG4Params;
			localG4Params.reserve(events.size());
			
			// Find the G4 hits on each surface
			for(const auto& event : events)
			{
				// Fast continue
				if(event.vertices().empty())
					continue;
				// Get the track ID that we follow
				const int trackID = event.vertices()[0]->particles_out()[0]->attribute<HepMC3::IntAttribute>("TrackID")->value();
				// The storage of each step
				std::vector<ActsExamples::SimParticle> g4Steps = collectG4Steps(event, trackID);
				
				const auto closestPoint = findClosestPoint(g4Steps, step.surface, gctx);
				if(closestPoint.has_value())
					localG4Params.push_back(closestPoint.value());
			}

			plot.scatter(localG4Params, surfaceParameters.eMeanPropagated, surfaceParameters.sMeanPropagated, surfaceParameters.surface);
			surfaceParameters.parametersG4 = localG4Params;
			// Calculate the mean of the G4 data
			const Acts::BoundVector meanG4 = calculateMean(localG4Params);
			surfaceParameters.meanG4 = std::move(meanG4);
			surfaceParameters.meanG4Free = Acts::detail::transformBoundToFreeParameters(*surfaceParameters.surface, gctx, surfaceParameters.meanG4);
			const Acts::BoundSymMatrix covG4 = calculateCovariance(localG4Params, meanG4);
			surfaceParameters.covG4 = std::move(covG4);
			
			summary.paramAtSurface.push_back(std::move(surfaceParameters));
		}
		summaries.push_back(std::move(summary));
	}
  // Plot the result
  plot.mean(summaries);

  return ActsExamples::ProcessCode::SUCCESS;
}

void
ActsExamples::MeanCalculator::plotMean(const std::vector<TrackSummary>& summaries) const {
	plot.mean(summaries);
	// TODO: should rather do this for all events once
	// TODO: split the plot in e.g. vs. r, phi, z, eta, ...
}
