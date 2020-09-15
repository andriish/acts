// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Fitter/KalmanFitterError.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <map>
#include <memory>

namespace Acts {

/// @brief Options struct how the Fitter is called
///
/// It contains the context of the fitter call, the outlier finder, the
/// optional surface where to express the fit result and configurations for
/// material effects and smoothing options
///
///
/// @note the context objects must be provided
template <typename outlier_finder_t = VoidOutlierFinder>
struct KalmanFitterOptions {
  // Broadcast the outlier finder type
  using OutlierFinder = outlier_finder_t;

  /// Deleted default constructor
  KalmanFitterOptions() = delete;

  /// PropagatorOptions with context
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param outlierFinder_ The outlier finder
  /// @param rSurface The reference surface for the fit to be expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param bwdFiltering Whether to run backward filtering as smoothing
  KalmanFitterOptions(std::reference_wrapper<const GeometryContext> gctx,
                      std::reference_wrapper<const MagneticFieldContext> mctx,
                      std::reference_wrapper<const CalibrationContext> cctx,
                      const OutlierFinder& outlierFinder_,
                      const Surface* rSurface = nullptr,
                      bool mScattering = true, bool eLoss = true,
                      bool bwdFiltering = false)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        outlierFinder(outlierFinder_),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        backwardFiltering(bwdFiltering) {}
    
  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The config for the outlier finder
  OutlierFinder outlierFinder;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering
  bool multipleScattering = true;

  /// Whether to consider energy loss
  bool energyLoss = true;

  /// Whether to run backward filtering
  bool backwardFiltering = false;
};

template <typename source_link_t>
struct KalmanFitterResult {
  // Fitted states that the actor has handled.
  MultiTrajectory<source_link_t> fittedStates;

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t trackTip = SIZE_MAX;

  // The optional Parameters at the provided surface
  std::optional<BoundParameters> fittedParameters;

  // Counter for states with measurements
  size_t measurementStates = 0;

  // Counter for handled states
  size_t processedStates = 0;

  // Indicator if smoothing has been done.
  bool smoothed = false;

  // Indicator if the propagation state has been reset
  bool reset = false;

  // Indicator if track fitting has been done
  bool finished = false;

  // Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  // Indicator if forward filtering has been done
  bool forwardFiltered = false;

  // Measurement surfaces handled in both forward and backward filtering
  std::vector<const GeometryObject*> passedAgainObject;

  // Iterator to the current volume and its source links
  const TrackingVolume* currentVolume = nullptr;
	
  Result<void> result{Result<void>::success()};
  
  struct FreeMeasurementCandidate
  {
	 const source_link_t* sourceLink = nullptr;
	 Vector3D position = Vector3D::Zero();
	 double distance = std::numeric_limits<double>::max();
  };
  std::vector<FreeMeasurementCandidate> currentFreeMeasurements;

};

/// @brief Kalman fitter implementation of Acts as a plugin
///
/// to the Propgator
///
/// @tparam propagator_t Type of the propagation class
/// @tparam updater_t Type of the kalman updater class
/// @tparam smoother_t Type of the kalman smoother class
/// @tparam outlier_finder_t Type of the outlier finder class
/// @tparam calibrator_t Type of the calibrator class
///
/// The Kalman filter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator
/// in order to initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update
/// and eventually the smoothing.  Updater, Smoother and Calibrator are
/// given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
/// - The Smoother is called at the end of the forward fit by the Actor.
/// - The outlier finder is called during the filtering by the Actor.
///   It determines if the measurement is an outlier
/// - The Calibrator is a dedicated calibration algorithm that allows
///   to calibrate measurements using track information, this could be
///    e.g. sagging for wires, module deformations, etc.
///
/// Measurements are not required to be ordered for the KalmanFilter,
/// measurement ordering needs to be figured out by the navigation of
/// the propagator.
///
/// The void components are provided mainly for unit testing.
template <typename propagator_t>
class KalmanFitter {
 public:
  /// The navigator type
  using KalmanNavigator = typename propagator_t::Navigator;
  using KalmanStepper = typename propagator_t::Stepper;
  
  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same<KalmanNavigator, DirectNavigator>::value;

  /// Default constructor is deleted
  KalmanFitter() = delete;

  /// Constructor from arguments
  KalmanFitter(propagator_t pPropagator,
               std::unique_ptr<const Logger> logger =
                   getDefaultLogger("KalmanFilter", Logging::INFO))
      : m_propagator(std::move(pPropagator)), m_logger(logger.release()) {}

 private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// Logger getter to support macros
  const Logger& logger() const { return *m_logger; }

  /// Owned logging instance
  std::shared_ptr<const Logger> m_logger;

 public:  
  /// @brief Propagator Actor plugin for the KalmanFilter
  ///
  /// @tparam source_link_t is an type fulfilling the @c SourceLinkConcept
  /// @tparam parameters_t The type of parameters used for "local" paremeters.
  ///
  /// The KalmanActor does not rely on the measurements to be
  /// sorted along the track.
  template <typename source_link_t, typename updater_t = VoidKalmanUpdater,
          typename smoother_t = VoidKalmanSmoother,
          typename outlier_finder_t = VoidOutlierFinder,
          typename calibrator_t = VoidMeasurementCalibrator>
  class Actor {
   public:
    using Self = Actor<source_link_t, updater_t, smoother_t, outlier_finder_t, calibrator_t>;
    using Fitter = KalmanFitter<propagator_t>;
    using Aborter = typename Fitter::template Aborter<Self>;
    using ActorList = typename std::conditional<isDirectNavigator, ActionList<DirectNavigator::Initializer, Self>, ActionList<Self>>::type;
    using Options = PropagatorOptions<ActorList, AbortList<Aborter, PathLimitReached>>;
    using State = typename propagator_t::template State<Options>;
	using SourceLink = source_link_t;
    using OutlierFinder = outlier_finder_t;
	
    /// Broadcast the result_type
    using result_type = KalmanFitterResult<source_link_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Allows retrieving measurements for a surface
    std::map<const GeometryObject*, source_link_t> boundInputMeasurements;
    /// Allows retrieving measurements for a volume
    std::map<const GeometryObject*, std::vector<source_link_t>> freeInputMeasurements;
	unsigned int numFreeInputMeasurements = 0;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether run smoothing as backward filtering
    bool backwardFiltering = false;

    /// @brief Kalman actor operation
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    void operator()(State& state, const KalmanStepper& stepper,
                    result_type& result) const {
      if (result.finished) {
        return;
      }

      ACTS_VERBOSE("KalmanFitter step");

      // This following is added due to the fact that the navigation
      // reinitialization in reverse call cannot guarantee the navigator to
      // target for extra layers in the backward-propagation starting volume.
      // Currently, manually set navigation stage to allow for targeting layers
      // after all the surfaces on the backward-propagation starting layer has
      // been processed. Otherwise, the navigation stage will be
      // Stage::boundaryTarget after navigator status call which means the extra
      // layers on the backward-propagation starting volume won't be targeted.
      // @Todo: Let the navigator do all the re-initialization
      if (result.reset and state.navigation.navSurfaceIter ==
                               state.navigation.navSurfaces.end()) {
        // So the navigator target call will target layers
        state.navigation.navigationStage = KalmanNavigator::Stage::layerTarget;
        // We only do this after the backward-propagation starting layer has
        // been processed
        result.reset = false;
      }

      // Update:
      // - Collect current free measurements
	  updateFreeMeasurementCandidates(state, stepper, result);
	  if(!result.currentFreeMeasurements.empty())
	  {
		  if(std::abs(result.currentFreeMeasurements[0].distance) < 1e-10)//state.stepping.tolerance)
		  {
			  // filter
			  if (state.stepping.navDir == forward and not result.smoothed and
				not result.forwardFiltered) {
				  filter(state, stepper, result);
				  result.currentFreeMeasurements.erase(result.currentFreeMeasurements.begin());
			  } 
			  else if (state.stepping.navDir == backward and
			   result.forwardFiltered) {
				 backwardFilter(state, stepper, result);
				 result.currentFreeMeasurements.erase(result.currentFreeMeasurements.begin());
			  }
			  if(!result.currentFreeMeasurements.empty())
			  {
				  calculateDistanceToMeasurement(state, stepper, result.currentFreeMeasurements[0]);
				  stepper.setStepSize(state.stepping, result.currentFreeMeasurements[0].distance * state.stepping.navDir, ConstrainedStep::user);
			  }
			  else
				  state.stepping.stepSize.release(ConstrainedStep::user);
		  }
		  else
		  {
			  stepper.setStepSize(state.stepping, result.currentFreeMeasurements[0].distance * state.stepping.navDir, ConstrainedStep::user);
		  }
	  }

	  // Update:
      // - Waiting for a current surface
      auto surface = state.navigation.currentSurface;
      if (surface != nullptr) {
        // Check if the surface is in the measurement map
        // -> Get the measurement / calibrate
        // -> Create the predicted state
        // -> Perform the kalman update
        // -> Check outlier behavior
        // -> Fill strack state information & update stepper information if
        // non-outlier
        if (state.stepping.navDir == forward and not result.smoothed and
            not result.forwardFiltered) {
          ACTS_VERBOSE("Perform forward filter step");
          auto res = filter(surface, state, stepper, result);
          if (!res.ok()) {
            ACTS_ERROR("Error in forward filter: " << res.error());
            result.result = res.error();
          }
        } else if (state.stepping.navDir == backward and
                   result.forwardFiltered) {
          ACTS_VERBOSE("Perform backward filter step");
          auto res = backwardFilter(surface, state, stepper, result);
          if (!res.ok()) {
            ACTS_ERROR("Error in backward filter: " << res.error());
            result.result = res.error();
          }
        }
      }

      // Finalization:
      // when all track states have been handled or the navigation is breaked,
      // reset navigation&stepping before run backward filtering or
      // proceed to run smoothing
      if (state.stepping.navDir == forward) {
        if ((result.measurementStates == boundInputMeasurements.size() + numFreeInputMeasurements) or
            (result.measurementStates > 0 and
             state.navigation.navigationBreak)) {
          if (backwardFiltering and not result.forwardFiltered) {
            ACTS_VERBOSE("Forward filtering done");
            result.forwardFiltered = true;
            // Start to run backward filtering:
            // Reverse navigation direction and reset navigation and stepping
            // state to last measurement
            ACTS_VERBOSE("Reverse navigation direction.");
            reverse(state, stepper, result);
          } else if (not result.smoothed) {
            // --> Search the starting state to run the smoothing
            // --> Call the smoothing
            // --> Set a stop condition when all track states have been
            // handled
            ACTS_VERBOSE("Finalize/run smoothing");
            auto res = finalize(state, stepper, result);
            if (!res.ok()) {
              ACTS_ERROR("Error in finalize: " << res.error());
              result.result = res.error();
            }
          }
        }
      }
      
      // Post-finalization:
      // - Progress to target/reference surface and built the final track
      // parameters
      if (result.smoothed or state.stepping.navDir == backward) {
        if (targetSurface == nullptr) {
          // If no target surface provided:
          // -> Return an error when using backward filtering mode
          // -> Fitting is finished here
          if (backwardFiltering) {
            ACTS_ERROR(
                "The target surface needed for aborting backward propagation "
                "is not provided");
            result.result =
                Result<void>(KalmanFitterError::BackwardUpdateFailed);
          } else {
            ACTS_VERBOSE(
                "No target surface set. Completing without fitted track "
                "parameter");
            // Remember the track fitting is done
            result.finished = true;
          }
        } else if (targetReached(state, stepper, *targetSurface)) {
          ACTS_VERBOSE("Completing with fitted track parameter");
          // Transport & bind the parameter to the final surface
          auto fittedState = stepper.boundState(state.stepping, *targetSurface);
          // Assign the fitted parameters
          result.fittedParameters = std::get<BoundParameters>(fittedState);

          // Reset smoothed status of states missed in backward filtering
          if (backwardFiltering) {
            result.fittedStates.applyBackwards(
                result.trackTip, [&](auto trackState) {
                  auto fSurface = &trackState.referenceObject();
                  auto surface_it = std::find_if(
                      result.passedAgainObject.begin(),
                      result.passedAgainObject.end(),
                      [=](const GeometryObject* s) { return s == fSurface; });
                  if (surface_it == result.passedAgainObject.end()) {
                    // If backward filtering missed this surface, then there is
                    // no smoothed parameter
                    trackState.data().iboundsmoothed =
                        detail_lt::IndexData::kInvalid;
                  }
                });
          }
          // Remember the track fitting is done
          result.finished = true;
        }
      }
    }

    /// @brief Kalman actor operation : reverse direction
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    void reverse(State& state, const KalmanStepper& stepper,
                 result_type& result) const {
      // Remember the navigation direciton has been reserved
      result.reset = true;

	// Reset FreeMeasurement components
	  result.currentVolume = nullptr;
	  result.currentFreeMeasurements.clear();
	  
      // Reset propagator options
      state.options.maxStepSize = -1.0 * state.options.maxStepSize;
      // Not sure if reset of pathLimit during propagation makes any sense
      state.options.pathLimit = -1.0 * state.options.pathLimit;

      // Reset stepping&navigation state using last measurement track state on
      // sensitive surface
      state.navigation = typename propagator_t::NavigatorState();
      result.fittedStates.applyBackwards(result.trackTip, [&](auto st) {
        if (st.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          const Surface* startSurface =
              dynamic_cast<const Surface*>(&st.referenceObject());
          if (startSurface != nullptr) {
            // Set the navigation state
            state.navigation.startSurface = startSurface;
            if (state.navigation.startSurface->associatedLayer() != nullptr) {
              state.navigation.startLayer =
                  state.navigation.startSurface->associatedLayer();
            }
            state.navigation.startVolume =
                state.navigation.startLayer->trackingVolume();
            state.navigation.targetSurface = targetSurface;
            state.navigation.currentSurface = state.navigation.startSurface;
            state.navigation.currentVolume = state.navigation.startVolume;

            // Update the stepping state
            stepper.resetState(state.stepping, st.boundFiltered(),
                               BoundMatrix(st.boundFilteredCovariance()), *startSurface, backward,
                               state.options.maxStepSize);

            // For the last measurement state, smoothed is filtered
            st.boundSmoothed() = st.boundFiltered();
            st.boundSmoothedCovariance() = st.boundFilteredCovariance();
            result.passedAgainObject.push_back(startSurface);

            // Update material effects for last measurement state in backward
            // direction
            materialInteractor(state.navigation.currentSurface, state, stepper);
          } else {
			const Volume* startVolume = dynamic_cast<const Volume*>(&st.referenceObject());
			// Set the navigation state
            state.navigation.startVolume = (TrackingVolume*) startVolume;
            state.navigation.startSurface = nullptr;
            state.navigation.startLayer =
				state.navigation.startVolume->associatedLayer(
                state.geoContext, stepper.position(state.stepping));
                
            state.navigation.targetSurface = targetSurface;
            state.navigation.currentSurface = state.navigation.startSurface;
            state.navigation.currentVolume = state.navigation.startVolume;
            
            // Update the stepping state
            stepper.resetState(state.stepping, st.freeFiltered(),
                               FreeMatrix(st.freeFilteredCovariance()), backward,
                               state.options.maxStepSize);

            // For the last measurement state, smoothed is filtered
            st.freeSmoothed() = st.freeFiltered();
            st.freeSmoothedCovariance() = st.freeFilteredCovariance();
            result.passedAgainObject.push_back(startVolume);
            
              updateFreeMeasurementCandidates(state, stepper, result);
			  if(!result.currentFreeMeasurements.empty())
			  {
				  if(std::abs(result.currentFreeMeasurements[0].distance) < state.stepping.tolerance)
				  {
					result.currentFreeMeasurements.erase(result.currentFreeMeasurements.begin());
					  if(!result.currentFreeMeasurements.empty())
					  {
						  calculateDistanceToMeasurement(state, stepper, result.currentFreeMeasurements[0]);
						  stepper.setStepSize(state.stepping, result.currentFreeMeasurements[0].distance * state.stepping.navDir, ConstrainedStep::user);
					  }
					  else
						  state.stepping.stepSize.release(ConstrainedStep::user);
				  }
				  else
				  {
					  // This shouldn't happen - the starting point should be within the tolerance
					  stepper.setStepSize(state.stepping, result.currentFreeMeasurements[0].distance * state.stepping.navDir, ConstrainedStep::user);
				  }
			  }
		  }
		  return false;  // abort execution
        }
        return true;  // continue execution
      });
    }

    /// @brief Kalman actor operation : update
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    Result<void> filter(const Surface* surface, State& state,
                        const KalmanStepper& stepper, result_type& result) const {
      // Try to find the surface in the measurement surfaces
      auto sourcelink_it = boundInputMeasurements.find(surface);
      if (sourcelink_it != boundInputMeasurements.end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geoID()
                                            << " detected.");

        // Transport & bind the state to the current surface
        auto [boundParams, jacobian, pathLength] =
            stepper.boundState(state.stepping, *surface);
            
        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, preUpdate);

		using TrackState = typename MultiTrajectory<source_link_t>::TrackStateProxy;
		TrackState trackStateProxy = std::visit([&](const auto& jac) -> TrackState {
			TrackStatePropMask mask = std::is_same<typename std::decay<decltype(jac)>::type, BoundMatrix>::value ? TrackStatePropMask::BoundAll : TrackStatePropMask::ToBoundAll;
        // add a full TrackState entry multi trajectory
        // (this allocates storage for all components, we will set them later)
        result.trackTip = result.fittedStates.addTrackState(
            mask, result.trackTip);

        // now get track state proxy back
        auto trackState =
            result.fittedStates.getTrackState(result.trackTip);
		// Store the jacobian
        if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, BoundMatrix>::value)
        {
        trackState.jacobianBoundToBound() = jac;
	}
	  if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, FreeToBoundMatrix>::value)
 {
		trackState.jacobianFreeToBound() = jac;
	}
      return trackState;}, jacobian);
	
        // assign the source link to the track state
        trackStateProxy.uncalibrated() = sourcelink_it->second;

        // Fill the track state
        trackStateProxy.boundPredicted() = boundParams.parameters();
        trackStateProxy.boundPredictedCovariance() = *boundParams.covariance();
        trackStateProxy.pathLength() = pathLength;
	
        // We have predicted parameters, so calibrate the uncalibrated input
        // measuerement
        std::visit(
            [&](const auto& calibrated) {
              trackStateProxy.setCalibrated(calibrated);
            },
            m_calibrator(trackStateProxy.uncalibrated(),
                         trackStateProxy.boundPredicted()));

        // Get and set the type flags
        auto& typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::MaterialFlag);
        typeFlags.set(TrackStateFlag::ParameterFlag);

        // Check if the state is an outlier.
        // If not, run Kalman update, tag it as a
        // measurement and update the stepping state. Otherwise, just tag it as
        // an outlier
        if (not m_outlierFinder(trackStateProxy)) {
          // Run Kalman update
          auto updateRes =
              m_updater(state.geoContext, trackStateProxy, forward);
          if (!updateRes.ok()) {
            ACTS_ERROR("Update step failed: " << updateRes.error());
            return updateRes.error();
          }
          // Set the measurement type flag
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE("Filtering step successful, updated parameters are : \n"
                       << trackStateProxy.boundFiltered().transpose());
          // update stepping state using filtered parameters after kalman
          stepper.update(state.stepping,
                         MultiTrajectoryHelpers::freeFiltered(
                             state.options.geoContext, trackStateProxy),
                         BoundSymMatrix(trackStateProxy.boundFilteredCovariance()));
          // We count the state with measurement
          ++result.measurementStates;
        } else {
          ACTS_VERBOSE(
              "Filtering step successful. But measurement is deterimined "
              "to "
              "be an outlier. Stepping state is not updated.")
          // Set the outlier type flag
          typeFlags.set(TrackStateFlag::OutlierFlag);
        }
        // We count the processed state
        ++result.processedStates;
      } else if (surface->surfaceMaterial() != nullptr) {
		  // TODO: there are probably some code duplications
        // We only create track states here if there is already measurement
        // detected
        if (result.measurementStates > 0) {
          if (surface->associatedDetectorElement() != nullptr) {
			// Transport & bind the state to the current surface
            auto [boundParams, jacobian, pathLength] =
                stepper.boundState(state.stepping, *surface);
			using TrackState = typename MultiTrajectory<source_link_t>::TrackStateProxy;
			TrackState trackStateProxy = std::visit([&](const auto& jac) -> TrackState {
				TrackStatePropMask mask = std::is_same<typename std::decay<decltype(jac)>::type, BoundMatrix>::value ? 
								TrackStatePropMask::BoundPredicted | TrackStatePropMask::BoundSmoothed | TrackStatePropMask::JacobianBoundToBound :
								TrackStatePropMask::BoundPredicted | TrackStatePropMask::BoundSmoothed | TrackStatePropMask::JacobianFreeToBound;
			// No source links on surface, add either hole or passive material
          // TrackState entry multi trajectory. No storage allocation for
          // uncalibrated/calibrated measurement and filtered parameter
			result.trackTip = result.fittedStates.addTrackState(
				mask, result.trackTip);

			// now get track state proxy back
			auto trackState =
				result.fittedStates.getTrackState(result.trackTip);
			// Store the jacobian
			if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, BoundMatrix>::value)
			{
				trackState.jacobianBoundToBound() = jac;
			}
			if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, FreeToBoundMatrix>::value)
			{
				trackState.jacobianFreeToBound() = jac;
			}
			return trackState;}, jacobian);

          // Set the surface
          trackStateProxy.setReferenceObject(surface->getSharedPtr());

          // Set the track state flags
          auto& typeFlags = trackStateProxy.typeFlags();
          typeFlags.set(TrackStateFlag::MaterialFlag);
          typeFlags.set(TrackStateFlag::ParameterFlag);
          
            ACTS_VERBOSE("Detected hole on " << surface->geoID());
            // If the surface is sensitive, set the hole type flag
            typeFlags.set(TrackStateFlag::HoleFlag);

            // Count the missed surface
            result.missedActiveSurfaces.push_back(surface);

            // Fill the track state
            trackStateProxy.boundPredicted() = boundParams.parameters();
            trackStateProxy.boundPredictedCovariance() =
                *boundParams.covariance();
            trackStateProxy.pathLength() = pathLength;
           
           // Set the filtered parameter index to be the same with predicted
          // parameter
          trackStateProxy.data().iboundfiltered =
              trackStateProxy.data().iboundpredicted;
              
          } else {
			// Transport & get curvilinear state instead of bound state
            auto [curvilinearParams, jacobian, pathLength] =
                stepper.curvilinearState(state.stepping);
			using TrackState = typename MultiTrajectory<source_link_t>::TrackStateProxy;
			TrackState trackStateProxy = std::visit([&](const auto& jac) -> TrackState {
				TrackStatePropMask mask = std::is_same<typename std::decay<decltype(jac)>::type, BoundMatrix>::value ? 
								TrackStatePropMask::BoundPredicted | TrackStatePropMask::BoundSmoothed | TrackStatePropMask::JacobianBoundToBound :
								TrackStatePropMask::BoundPredicted | TrackStatePropMask::BoundSmoothed | TrackStatePropMask::JacobianFreeToBound;
			// No source links on surface, add either hole or passive material
          // TrackState entry multi trajectory. No storage allocation for
          // uncalibrated/calibrated measurement and filtered parameter
			result.trackTip = result.fittedStates.addTrackState(
				mask, result.trackTip);

			// now get track state proxy back
			auto trackState =
				result.fittedStates.getTrackState(result.trackTip);
			// Store the jacobian
			if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, BoundMatrix>::value)
			{
				trackState.jacobianBoundToBound() = jac;
			}
			if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, FreeToBoundMatrix>::value)
			{
				trackState.jacobianFreeToBound() = jac;
			}
			return trackState;}, jacobian);

          // Set the surface
          trackStateProxy.setReferenceObject(surface->getSharedPtr());

          // Set the track state flags
          auto& typeFlags = trackStateProxy.typeFlags();
          typeFlags.set(TrackStateFlag::MaterialFlag);
          typeFlags.set(TrackStateFlag::ParameterFlag);
          
            ACTS_VERBOSE("Detected in-sensitive surface " << surface->geoID());

            // Fill the track state
            trackStateProxy.boundPredicted() = curvilinearParams.parameters();
            trackStateProxy.boundPredictedCovariance() =
                *curvilinearParams.covariance();
            trackStateProxy.pathLength() = pathLength;
            
           // Set the filtered parameter index to be the same with predicted
          // parameter
          trackStateProxy.data().iboundfiltered =
              trackStateProxy.data().iboundpredicted;
          }

          // We count the processed state
          ++result.processedStates;
        }

        // Update state and stepper with material effects
        materialInteractor(surface, state, stepper, fullUpdate);
      }
      return Result<void>::success();
    }

	/// @brief Kalman actor operation : update
    ///
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    Result<void> filter(State& state, const KalmanStepper& stepper, result_type& result) const {
	  const source_link_t& sourceLink = *result.currentFreeMeasurements[0].sourceLink;

		auto [freeParams, jacobian, pathLength] = stepper.freeState(state.stepping);

		using TrackState = typename MultiTrajectory<source_link_t>::TrackStateProxy;
		TrackState trackStateProxy = std::visit([&](const auto& jac) -> TrackState {
			TrackStatePropMask mask = std::is_same<typename std::decay<decltype(jac)>::type, FreeMatrix>::value ? TrackStatePropMask::FreeAll : TrackStatePropMask::ToFreeAll;
        // add a full TrackState entry multi trajectory
        // (this allocates storage for all components, we will set them later)
        result.trackTip = result.fittedStates.addTrackState(
            mask, result.trackTip);

        // now get track state proxy back
        auto trackState =
            result.fittedStates.getTrackState(result.trackTip);

		// Store the jacobian
        if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, FreeMatrix>::value)
        {
        trackState.jacobianFreeToFree() = jac;
	}
	  if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, BoundToFreeMatrix>::value)
 {
		trackState.jacobianBoundToFree() = jac;
	}
      return trackState;}, jacobian);

        // assign the source link to the track state
        trackStateProxy.uncalibrated() = sourceLink;

        // Fill the track state
        trackStateProxy.freePredicted() = freeParams.parameters();
        trackStateProxy.freePredictedCovariance() = *freeParams.covariance();
        trackStateProxy.pathLength() = pathLength;

        // We have predicted parameters, so calibrate the uncalibrated input
        // measuerement
        std::visit(
            [&](const auto& calibrated) {
              trackStateProxy.setCalibrated(calibrated);
            },
            m_calibrator(trackStateProxy.uncalibrated(),
                         trackStateProxy.freePredicted()));

        // Get and set the type flags
        auto& typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::MaterialFlag);
        typeFlags.set(TrackStateFlag::ParameterFlag);
        // If the update is successful, set covariance and
        auto updateRes = m_updater(state.geoContext, trackStateProxy, forward, false);
        if (!updateRes.ok()) {
          ACTS_ERROR("Update step failed: " << updateRes.error());
          return updateRes.error();
        } else {
          if (not m_outlierFinder(trackStateProxy)) {
            // Set the measurement type flag
            typeFlags.set(TrackStateFlag::MeasurementFlag);
            // Update the stepping state with filtered parameters
            ACTS_VERBOSE(
                "Filtering step successful, updated parameters are : \n"
                << trackStateProxy.freeFiltered().transpose());
            // update stepping state using filtered parameters after kalman
            // update We need to (re-)construct a BoundParameters instance
            // here, which is a bit awkward.
            stepper.update(
                state.stepping,
                trackStateProxy.freeFiltered(),
                FreeSymMatrix(trackStateProxy.freeFilteredCovariance()));
            // We count the state with measurement
            ++result.measurementStates;
          } else {
            ACTS_VERBOSE(
                "Filtering step successful. But measurement is deterimined "
                "to "
                "be an outlier. Stepping state is not updated.")
            // Set the outlier type flag
            typeFlags.set(TrackStateFlag::OutlierFlag);
          }
        }
        // We count the processed state
        ++result.processedStates;
      return Result<void>::success();
    }
    
    /// @brief Kalman actor operation : backward update
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    Result<void> backwardFilter(const Surface* surface,
                                State& state,
                                const KalmanStepper& stepper,
                                result_type& result) const {
      // Try to find the surface in the measurement surfaces
      auto sourcelink_it = boundInputMeasurements.find(surface);
      if (sourcelink_it != boundInputMeasurements.end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface "
                     << surface->geoID()
                     << " detected in backward propagation.");

        // No backward filtering for last measurement state, but still update
        // with material effects
        if (state.stepping.navDir == backward and
            surface == state.navigation.startSurface) {
          materialInteractor(surface, state, stepper);
          return Result<void>::success();
        }

        // Transport & bind the state to the current surface
        auto [boundParams, jacobian, pathLength] =
            stepper.boundState(state.stepping, *surface);
            
        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, preUpdate);

		using TrackState = typename MultiTrajectory<source_link_t>::TrackStateProxy;
		TrackState trackStateProxy = std::visit([&](const auto& jac) -> TrackState {
			TrackStatePropMask mask = std::is_same<typename std::decay<decltype(jac)>::type, BoundMatrix>::value ? TrackStatePropMask::BoundAll : TrackStatePropMask::ToBoundAll;
        // Create a detached track state proxy
        auto tempTrackTip = result.fittedStates.addTrackState(mask);

        // now get track state proxy back
        auto trackState =
            result.fittedStates.getTrackState(tempTrackTip);
		// Store the jacobian
        if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, BoundMatrix>::value)
        {
			trackState.jacobianBoundToBound() = jac;
		}
		if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, FreeToBoundMatrix>::value)
		{
			trackState.jacobianFreeToBound() = jac;
		}
		return trackState;}, jacobian);

        // Assign the source link to the detached track state
        trackStateProxy.uncalibrated() = sourcelink_it->second;

        // Fill the track state
        trackStateProxy.boundPredicted() = boundParams.parameters();
        trackStateProxy.boundPredictedCovariance() = *boundParams.covariance();
        trackStateProxy.pathLength() = pathLength;

        // We have predicted parameters, so calibrate the uncalibrated input
        // measuerement
        std::visit(
            [&](const auto& calibrated) {
              trackStateProxy.setCalibrated(calibrated);
            },
            m_calibrator(trackStateProxy.uncalibrated(),
                         trackStateProxy.boundPredicted()));
        // If the update is successful, set covariance and
        auto updateRes = m_updater(state.geoContext, trackStateProxy, backward);
        if (!updateRes.ok()) {
          ACTS_ERROR("Backward update step failed: " << updateRes.error());
          return updateRes.error();
        } else {
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE(
              "Backward Filtering step successful, updated parameters are : "
              "\n"
              << trackStateProxy.boundFiltered().transpose());

          // Fill the smoothed parameter for the existing track state
          result.fittedStates.applyBackwards(
              result.trackTip, [&](auto trackState) {
                auto fReferenceObject = &trackState.referenceObject();
                if (fReferenceObject == surface) {
                  result.passedAgainObject.push_back(surface);
                  trackState.boundSmoothed() = trackStateProxy.boundFiltered();
                  trackState.boundSmoothedCovariance() =
                      trackStateProxy.boundFilteredCovariance();
                  return false;
                }
                return true;
              });

          // update stepping state using filtered parameters after kalman
          // update
          stepper.update(state.stepping,
                         MultiTrajectoryHelpers::freeFiltered(
                             state.options.geoContext, trackStateProxy),
                         BoundSymMatrix(trackStateProxy.boundFilteredCovariance()));

          // Update state and stepper with post material effects
          materialInteractor(surface, state, stepper, postUpdate);
        }
      } else if (surface->surfaceMaterial() != nullptr) {
        // Transport covariance
        if (surface->associatedDetectorElement() != nullptr) {
          ACTS_VERBOSE("Detected hole on " << surface->geoID()
                                           << " in backward filtering");
          if (state.stepping.covTransport) {
            stepper.covarianceTransport(state.stepping, *surface);
          }
        } else {
          ACTS_VERBOSE("Detected in-sensitive surface "
                       << surface->geoID() << " in backward filtering");
          if (state.stepping.covTransport) {
            stepper.covarianceTransport(state.stepping);
          }
        }

        // Update state and stepper with material effects
        materialInteractor(surface, state, stepper);
      }
      return Result<void>::success();
    }
    
    /// @brief Kalman actor operation : backward update
    ///
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    Result<void> backwardFilter(State& state,
                                const KalmanStepper& stepper,
                                result_type& result) const {
	const source_link_t& sourceLink = *result.currentFreeMeasurements[0].sourceLink;

		auto [freeParams, jacobian, pathLength] = stepper.freeState(state.stepping);

		using TrackState = typename MultiTrajectory<source_link_t>::TrackStateProxy;
		TrackState trackStateProxy = std::visit([&](const auto& jac) -> TrackState {
			TrackStatePropMask mask = std::is_same<typename std::decay<decltype(jac)>::type, FreeMatrix>::value ? TrackStatePropMask::FreeAll : TrackStatePropMask::ToFreeAll;
        // Create a detached track state proxy
        auto tempTrackTip = result.fittedStates.addTrackState(mask);

        // now get track state proxy back
        auto trackState = result.fittedStates.getTrackState(tempTrackTip);
		// Store the jacobian
        if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, FreeMatrix>::value)
        {
			trackState.jacobianFreeToFree() = jac;
		}
		if constexpr (std::is_same<typename std::decay<decltype(jac)>::type, BoundToFreeMatrix>::value)
		{
			trackState.jacobianBoundToFree() = jac;
		}
		return trackState;}, jacobian);

        // Assign the source link to the detached track state
        trackStateProxy.uncalibrated() = sourceLink;

        // Fill the track state
        trackStateProxy.freePredicted() = freeParams.parameters();
        trackStateProxy.freePredictedCovariance() = *freeParams.covariance();
        trackStateProxy.pathLength() = pathLength;

        // We have predicted parameters, so calibrate the uncalibrated input
        // measuerement
        std::visit(
            [&](const auto& calibrated) {
              trackStateProxy.setCalibrated(calibrated);
            },
            m_calibrator(trackStateProxy.uncalibrated(),
                         trackStateProxy.freePredicted()));

        // If the update is successful, set covariance and
        auto updateRes = m_updater(state.geoContext, trackStateProxy, backward, false);
        if (!updateRes.ok()) {
          ACTS_ERROR("Backward update step failed: " << updateRes.error());
          return updateRes.error();
        } else {
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE(
              "Backward Filtering step successful, updated parameters are : "
              "\n"
              << trackStateProxy.freeFiltered().transpose());

          // Fill the smoothed parameter for the existing track state
          result.fittedStates.applyBackwards(
              result.trackTip, [&](auto trackState) {
                auto fReferenceObject = &trackState.referenceObject();
                if (trackState.uncalibrated() == trackStateProxy.uncalibrated()) {
                  result.passedAgainObject.push_back(fReferenceObject);
                  trackState.freeSmoothed() = trackStateProxy.freeFiltered();
                  trackState.freeSmoothedCovariance() =
                      trackStateProxy.freeFilteredCovariance();
                  return false;
                }
                return true;
              });

          // update stepping state using filtered parameters after kalman
          // update
            stepper.update(
                state.stepping,
                trackStateProxy.freeFiltered(),
                FreeSymMatrix(trackStateProxy.freeFilteredCovariance()));
        }
      return Result<void>::success();
    }

    /// @brief Kalman actor operation : material interaction
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param updateStage The materal update stage
    ///
    void materialInteractor(
        const Surface* surface, State& state, const KalmanStepper& stepper,
        const MaterialUpdateStage& updateStage = fullUpdate) const {
      // Indicator if having material
      bool hasMaterial = false;

      if (surface and surface->surfaceMaterial()) {
        // Prepare relevant input particle properties
        detail::PointwiseMaterialInteraction interaction(surface, state,
                                                         stepper);
        // Evaluate the material properties
        if (interaction.evaluateMaterialProperties(state, updateStage)) {
          // Surface has material at this stage
          hasMaterial = true;

          // Evaluate the material effects
          interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                           energyLoss);

          // Screen out material effects info
          ACTS_VERBOSE("Material effects on surface: "
                       << surface->geoID()
                       << " at update stage: " << updateStage << " are :");
          ACTS_VERBOSE("eLoss = "
                       << interaction.Eloss << ", "
                       << "variancePhi = " << interaction.variancePhi << ", "
                       << "varianceTheta = " << interaction.varianceTheta
                       << ", "
                       << "varianceQoverP = " << interaction.varianceQoverP);
          // Update the state and stepper with material effects
          interaction.updateState(state, stepper);
        }
      }

      if (not hasMaterial) {
        // Screen out message
        ACTS_VERBOSE("No material effects on surface: " << surface->geoID()
                                                        << " at update stage: "
                                                        << updateStage);
      }
    }

    /// @brief Kalman actor operation : finalize
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    Result<void> finalize(State& state, const KalmanStepper& stepper,
                          result_type& result) const {
      // Remember you smoothed the track states
      result.smoothed = true;

      // Get the indices of measurement states;
      std::vector<size_t> measurementIndices;
      measurementIndices.reserve(result.measurementStates);
      // Count track states to be smoothed
      size_t nStates = 0;
      result.fittedStates.applyBackwards(result.trackTip, [&](auto st) {
        bool isMeasurement =
            st.typeFlags().test(TrackStateFlag::MeasurementFlag);
        if (isMeasurement) {
          measurementIndices.emplace_back(st.index());
        } else if (measurementIndices.empty()) {
          // No smoothed parameters if the last measurement state has not been
          // found yet
          st.data().iboundsmoothed = detail_lt::IndexData::kInvalid;
        }
        // Start count when the last measurement state is found
        if (not measurementIndices.empty()) {
          nStates++;
        }
      });
      // Return error if the track has no measurement states (but this should
      // not happen)
      if (measurementIndices.empty()) {
        ACTS_ERROR("Smoothing for a track without measurements.");
        return KalmanFitterError::SmoothFailed;
      }
      // Screen output for debugging
      if (logger().doPrint(Logging::VERBOSE)) {
        ACTS_VERBOSE("Apply smoothing on " << nStates
                                           << " filtered track states.");
      }

      // Smooth the track states
      auto smoothRes = m_smoother(state.geoContext, result.fittedStates,
                                  measurementIndices.front());
      if (!smoothRes.ok()) {
        ACTS_ERROR("Smoothing step failed: " << smoothRes.error());
        return smoothRes.error();
      }
      // Obtain the smoothed parameters at first measurement state
      auto firstMeasurement =
          result.fittedStates.getTrackState(measurementIndices.back());

      // Update the stepping parameters - in order to progress to destination
      ACTS_VERBOSE(
          "Smoothing successful, updating stepping state, "
          "set target surface.");
          if(firstMeasurement.hasBoundSmoothed())
          {
			stepper.update(state.stepping,
                     MultiTrajectoryHelpers::freeSmoothed(
                         state.options.geoContext, firstMeasurement),
                     BoundSymMatrix(firstMeasurement.boundSmoothedCovariance()));
           }   
           else
			   stepper.update(state.stepping,
							 firstMeasurement.freeSmoothed(),
						 FreeSymMatrix(firstMeasurement.freeSmoothedCovariance()));
      // Reverse the propagation direction
      state.stepping.stepSize =
          ConstrainedStep(-1. * state.options.maxStepSize);
      state.stepping.navDir = backward;
      // Set accumulatd path to zero before targeting surface
      state.stepping.pathAccumulated = 0.;
      return Result<void>::success();
    }

	void updateFreeMeasurementCandidates(State& state, const KalmanStepper& stepper, result_type& result) const
	{
		// Test whether the volume changed
      if(result.currentVolume == nullptr || state.navigation.currentVolume != result.currentVolume)
      {
		  result.currentVolume = state.navigation.currentVolume;
		  // Search for measurements in the volume
		  auto measurementsInCurrentVolume = freeInputMeasurements.find(state.navigation.currentVolume);
		  if(measurementsInCurrentVolume != freeInputMeasurements.end())
		  {
			 result.currentFreeMeasurements = collectCandidates(measurementsInCurrentVolume->second);
			 calculateAllDistancesToMeasurement(state, stepper, result.currentFreeMeasurements);
			 sortFreeInputMeasurements(result.currentFreeMeasurements);
			 while(!result.currentFreeMeasurements.empty() && result.currentFreeMeasurements[0].distance < 0.)
			 {
				result.currentFreeMeasurements.erase(result.currentFreeMeasurements.begin()); 
			 }
		  }
		  else
		  {
			  state.stepping.stepSize.release(ConstrainedStep::user);
			  result.currentFreeMeasurements.clear();
		  }
	  }
	  else
	  {
		  // Update the distance to the closest element
		  if(!result.currentFreeMeasurements.empty())
			calculateDistanceToMeasurement(state, stepper, result.currentFreeMeasurements[0]);
	  }
	 }
	  

  std::vector<typename result_type::FreeMeasurementCandidate> 
  sortFreeInputMeasurements(std::vector<typename result_type::FreeMeasurementCandidate>& candidates) const {	
	// Sort elements by distance
	std::sort(candidates.begin(), candidates.end(), [&](const typename result_type::FreeMeasurementCandidate& a, const typename result_type::FreeMeasurementCandidate& b) 
		{ return a.distance < b.distance;});
    return candidates;
  }
  
  std::vector<typename result_type::FreeMeasurementCandidate> 
  collectCandidates(const std::vector<source_link_t>& sourceLinks) const {
	  const unsigned int numberOfElements = sourceLinks.size();
	  std::vector<typename result_type::FreeMeasurementCandidate> candidates(numberOfElements);

	for(unsigned int i = 0; i < numberOfElements; i++)
	{
		std::visit([&](const auto& meas) {
		   using MeasType = typename std::decay<decltype(meas)>::type;
		   using Indices = typename MeasType::ParametersIndices;
		   if constexpr (std::is_same<Indices, FreeParametersIndices>::value)
		   {
				// Store the position
				candidates[i].position = (meas.projector().transpose() * meas.parameters()).template segment<3>(eFreePos0);
			  
				// Store source link
				candidates[i].sourceLink = &sourceLinks[i];
			}
		}, *sourceLinks[i]);
	}
	return candidates;
  }
  
  void
  calculateAllDistancesToMeasurement(const State& state, const KalmanStepper& stepper, std::vector<typename result_type::FreeMeasurementCandidate>& candidates) const
  {
	// Evaluate all distances	     	  
	for(unsigned int i = 0; i < candidates.size(); i++)
	{
		calculateDistanceToMeasurement(state, stepper, candidates[i]);
	} 
  }
  
  void
  calculateDistanceToMeasurement(const State& state, const KalmanStepper& stepper, typename result_type::FreeMeasurementCandidate& candidate) const
  {
	using Plane = Eigen::Hyperplane<double, 3>;
	
	 // Get the position & direction
    auto position = stepper.position(state.stepping);
    auto direction = stepper.direction(state.stepping);
	
	// Build the plane
	auto plane = Plane(direction, candidate.position);
  
	// Store intersection distance and source link
	candidate.distance = -1. * state.stepping.navDir * plane.signedDistance(position);
	candidate.distance = candidate.position.x() - position.x();
  }
  
    /// Pointer to a logger that is owned by the parent, KalmanFilter
    const Logger* m_logger;

    /// Getter for the logger, to support logging macros
    const Logger& logger() const { return *m_logger; }

    /// The Kalman updater
    updater_t m_updater;

    /// The Kalman smoother
    smoother_t m_smoother;

    /// The outlier finder
    outlier_finder_t m_outlierFinder;

    /// The Measuremetn calibrator
    calibrator_t m_calibrator;

    /// The Surface beeing
    SurfaceReached targetReached;
  };

  template <typename actor_type_t>
  class Aborter {
   public:
    /// Broadcast the result_type
    using action_type = actor_type_t;
    
    bool operator()(typename action_type::State& /*state*/, const KalmanStepper& /*stepper*/,
                    const typename action_type::result_type& result) const {
      if (!result.result.ok() or result.finished) {
        return true;
      }
      return false;
    }
  };

private:
    template <typename actor_type_t>
    auto
    buildPropagatorOptions(const std::vector<typename actor_type_t::SourceLink>& boundSourcelinks, const std::vector<typename actor_type_t::SourceLink>& freeSourcelinks, 
    const KalmanFitterOptions<typename actor_type_t::OutlierFinder>& kfOptions) const
    {
		using SourceLink = typename actor_type_t::SourceLink;
        static_assert(SourceLinkConcept<SourceLink>,
                  "Source link does not fulfill SourceLinkConcept");

		// To be able to find measurements later, we put them into a map
		// We need to copy input SourceLinks anyways, so the map can own them.
		ACTS_VERBOSE("Preparing " << boundSourcelinks.size() << " bound input measurements");
		std::map<const GeometryObject*, SourceLink> boundInputMeasurements;
		for (const auto& sl : boundSourcelinks) {
		  const GeometryObject* srf = &sl.referenceObject();
		  boundInputMeasurements.emplace(srf, sl);
		}
		ACTS_VERBOSE("Preparing " << freeSourcelinks.size() << " free input measurements");
		std::map<const GeometryObject*, std::vector<SourceLink>> freeInputMeasurements;
		for (const auto& sl : freeSourcelinks) {
		  const GeometryObject* srf = &sl.referenceObject();
		  freeInputMeasurements[srf].emplace_back(sl);
		}

		// Create the ActionList and AbortList
		using KalmanAborter = Aborter<actor_type_t>;
		using KalmanActor = actor_type_t;
		using Aborters = AbortList<KalmanAborter>;
		using Actors = typename std::conditional<isDirectNavigator, ActionList<DirectNavigator::Initializer, KalmanActor>, ActionList<KalmanActor>>::type;
		// Create relevant options for the propagation options
		PropagatorOptions<Actors, Aborters> kalmanOptions(
			kfOptions.geoContext, kfOptions.magFieldContext);

		// Catch the actor and set the measurements
		auto& kalmanActor = kalmanOptions.actionList.template get<KalmanActor>();
		kalmanActor.m_logger = m_logger.get();
		kalmanActor.boundInputMeasurements = std::move(boundInputMeasurements);
		kalmanActor.freeInputMeasurements = std::move(freeInputMeasurements);
		kalmanActor.numFreeInputMeasurements = freeSourcelinks.size();
		kalmanActor.targetSurface = kfOptions.referenceSurface;
		kalmanActor.multipleScattering = kfOptions.multipleScattering;
		kalmanActor.energyLoss = kfOptions.energyLoss;
		kalmanActor.backwardFiltering = kfOptions.backwardFiltering;

		// Set config for outlier finder
		kalmanActor.m_outlierFinder = kfOptions.outlierFinder;

		// also set logger on updater and smoother
		kalmanActor.m_updater.m_logger = m_logger;
		kalmanActor.m_smoother.m_logger = m_logger;
				
		return kalmanOptions;
	}

public:   
  /// Fit implementation of the foward filter, calls the
  /// the forward filter and backward smoother
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @note The input measurements are given in the form of @c SourceLinks.
  /// It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename updater_t = VoidKalmanUpdater,
          typename smoother_t = VoidKalmanSmoother,
          typename calibrator_t = VoidMeasurementCalibrator, typename outlier_finder_t, typename source_link_t, typename start_parameters_t>
  auto fit(const std::vector<source_link_t>& sourcelinks,
           const start_parameters_t& sParameters,
           const KalmanFitterOptions<outlier_finder_t>& kfOptions,
           const std::vector<source_link_t>& freeSourcelinks = {}) const
      -> std::enable_if_t<!isDirectNavigator,
                          Result<KalmanFitterResult<source_link_t>>> {
	using ActorType = Actor<source_link_t, updater_t, smoother_t, outlier_finder_t, calibrator_t>;
    using KalmanResult = typename ActorType::result_type;
   
    // Create relevant options for the propagation options
    auto kalmanOptions = buildPropagatorOptions<ActorType>(sourcelinks, freeSourcelinks, kfOptions);

    // Run the fitter
    auto result = m_propagator.template propagate(sParameters, kalmanOptions);

    if (!result.ok()) {
      ACTS_ERROR("Propapation failed: " << result.error());
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the fit
    auto kalmanResult = propRes.template get<KalmanResult>();

    /// It could happen that the fit ends in zero processed states.
    /// The result gets meaningless so such case is regarded as fit failure.
    if (kalmanResult.result.ok() and not kalmanResult.processedStates) {
      kalmanResult.result = Result<void>(KalmanFitterError::NoMeasurementFound);
    }

    if (!kalmanResult.result.ok()) {
      ACTS_ERROR("KalmanFilter failed: " << kalmanResult.result.error());
      return kalmanResult.result.error();
    }

    // Return the converted Track
    return kalmanResult;
  }

  /// Fit implementation of the foward filter, calls the
  /// the forward filter and backward smoother
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @param sSequence surface sequence used to initialize a DirectNavigator
  /// @note The input measurements are given in the form of @c SourceLinks.
  /// It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename updater_t = VoidKalmanUpdater,
          typename smoother_t = VoidKalmanSmoother,
          typename calibrator_t = VoidMeasurementCalibrator, typename outlier_finder_t, typename source_link_t, typename start_parameters_t>
  auto fit(const std::vector<source_link_t>& sourcelinks,
           const start_parameters_t& sParameters,
           const KalmanFitterOptions<outlier_finder_t>& kfOptions,
           const std::vector<const Surface*>& sSequence) const
      -> std::enable_if_t<isDirectNavigator,
                          Result<KalmanFitterResult<source_link_t>>> {
	using ActorType = Actor<source_link_t, updater_t, smoother_t, outlier_finder_t, calibrator_t>;
    using KalmanResult = typename ActorType::result_type;
   
    // Create relevant options for the propagation options
    auto kalmanOptions = buildPropagatorOptions<ActorType>(sourcelinks, {}, kfOptions);

    // Set the surface sequence
    auto& dInitializer =
        kalmanOptions.actionList.template get<DirectNavigator::Initializer>();
    dInitializer.surfaceSequence = sSequence;

    // Run the fitter
    auto result = m_propagator.template propagate(sParameters, kalmanOptions);

    if (!result.ok()) {
      ACTS_ERROR("Propapation failed: " << result.error());
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the fit
    auto kalmanResult = propRes.template get<KalmanResult>();

    /// It could happen that the fit ends in zero processed states.
    /// The result gets meaningless so such case is regarded as fit failure.
    if (kalmanResult.result.ok() and not kalmanResult.processedStates) {
      kalmanResult.result = Result<void>(KalmanFitterError::NoMeasurementFound);
    }

    if (!kalmanResult.result.ok()) {
      ACTS_ERROR("KalmanFilter failed: " << kalmanResult.result.error());
      return kalmanResult.result.error();
    }

    // Return the converted Track
    return kalmanResult;
  }
};

}  // namespace Acts
