// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/covariance_helper.hpp"
#include "Acts/Fitter/KalmanFitterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>

#include <boost/range/adaptors.hpp>

namespace Acts {

/// @brief Kalman smoother implementation based on Gain matrix formalism
///
/// @tparam parameters_t Type of the track parameters
/// @tparam jacobian_t Type of the Jacobian
class GainMatrixSmoother {
 public:
  /// @brief Gain Matrix smoother implementation
  ///

  /// Constructor with (non-owning) logger
  /// @param logger a logger instance
  GainMatrixSmoother(
      std::shared_ptr<const Logger> logger = std::shared_ptr<const Logger>(
          getDefaultLogger("GainMatrixSmoother", Logging::INFO).release()));

  /// Operater for Kalman smoothing
  ///
  /// @tparam source_link_t The type of source link
  ///
  /// @param gctx The geometry context for the smoothing
  /// @param trajectory The trajectory to be smoothed
  /// @param entryIndex The index of state to start the smoothing
  /// @param globalTrackParamsCovPtr The pointer to global track parameters
  /// covariance matrix
  ///
  /// @return The smoothed track parameters at the first measurement state
  template <typename source_link_t>
  Result<void> operator()(const GeometryContext& /* gctx */,
                          MultiTrajectory<source_link_t>& trajectory,
                          size_t entryIndex) const {
	ACTS_VERBOSE("Invoked GainMatrixSmoother on entry index: " << entryIndex);
	using namespace boost::adaptors;

	// For the last state: smoothed is filtered - also: switch to next
	ACTS_VERBOSE("Getting previous track state");
	auto prev_ts = trajectory.getTrackState(entryIndex);

	if(prev_ts.hasBoundFiltered())
	{
		prev_ts.boundSmoothed() = prev_ts.boundFiltered();
		prev_ts.boundSmoothedCovariance() = prev_ts.boundFilteredCovariance();
	}
	else
	{
		prev_ts.freeSmoothed() = prev_ts.freeFiltered();
		prev_ts.freeSmoothedCovariance() = prev_ts.freeFilteredCovariance();  				
	}

	// make sure there is more than one track state
	std::optional<std::error_code> error{std::nullopt};  // assume ok
	if (prev_ts.previous() == Acts::detail_lt::IndexData::kInvalid) {
	  ACTS_VERBOSE("Only one track state given, smoothing terminates early");
	} else {
	  ACTS_VERBOSE("Start smoothing from previous track state at index: "
				   << prev_ts.previous());

	  trajectory.applyBackwards(prev_ts.previous(), [&prev_ts, &error,
													 this](auto ts) {
		// should have filtered and predicted, this should also include the
		// covariances.
		if(ts.hasBoundPredicted())
		{
			assert(ts.hasBoundFiltered());
			// previous trackstate should have smoothed and predicted
			if(prev_ts.hasBoundPredicted())
			{
				assert(prev_ts.hasBoundSmoothed());
				assert(prev_ts.hasJacobianBoundToBound());
				if(smooth<eBoundParametersSize, eBoundParametersSize>(error, ts.boundFiltered(), ts.boundSmoothed(), ts.boundFilteredCovariance(), ts.boundSmoothedCovariance(),
																				prev_ts.boundPredicted(), prev_ts.boundSmoothed(), prev_ts.jacobianBoundToBound(),
																				prev_ts.boundPredictedCovariance(), prev_ts.boundSmoothedCovariance()))
				{
					prev_ts = ts;
					return true;  // continue execution
				}
				else
					return false;		
			}
			else
			{
				assert(prev_ts.hasFreeSmoothed());
				assert(prev_ts.hasFreePredicted());
				assert(prev_ts.hasJacobianBoundToFree());
				if(smooth<eFreeParametersSize, eBoundParametersSize>(error, ts.boundFiltered(), ts.boundSmoothed(), ts.boundFilteredCovariance(), ts.boundSmoothedCovariance(),
																				prev_ts.freePredicted(), prev_ts.freeSmoothed(), prev_ts.jacobianBoundToFree(),
																				prev_ts.freePredictedCovariance(), prev_ts.freeSmoothedCovariance()))
				{
					prev_ts = ts;
					return true;  // continue execution
				}
				else
					return false;
			}	
		}
		else
		{
			assert(ts.hasFreeFiltered());
			assert(ts.hasFreePredicted());
			assert(prev_ts.hasJacobianFreeToFree() || prev_ts.hasJacobianBoundToFree());

			// previous trackstate should have smoothed and predicted
			if(prev_ts.hasBoundPredicted())
			{
				assert(prev_ts.hasBoundSmoothed());
				assert(prev_ts.hasJacobianFreeToBound());
				if(smooth<eBoundParametersSize, eFreeParametersSize>(error, ts.freeFiltered(), ts.freeSmoothed(), ts.freeFilteredCovariance(), ts.freeSmoothedCovariance(),
																				prev_ts.boundPredicted(), prev_ts.boundSmoothed(), prev_ts.jacobianFreeToBound(),
																				prev_ts.boundPredictedCovariance(), prev_ts.boundSmoothedCovariance()))
				{
					prev_ts = ts;
					return true;  // continue execution
				}
				else
					return false;
			}
			else
			{
				assert(prev_ts.hasFreeSmoothed());
				assert(prev_ts.hasFreePredicted());
				assert(prev_ts.hasJacobianFreeToFree());
				if(smooth<eFreeParametersSize, eFreeParametersSize>(error, ts.freeFiltered(), ts.freeSmoothed(), ts.freeFilteredCovariance(), ts.freeSmoothedCovariance(),
																				prev_ts.freePredicted(), prev_ts.freeSmoothed(), prev_ts.jacobianFreeToFree(),
																				prev_ts.freePredictedCovariance(), prev_ts.freeSmoothedCovariance()))
				{
					prev_ts = ts;
					return true;  // continue execution
				}
				else
					return false;
			}

		}
	  });
	}
    if (error) {
      // error is set, return result
      return *error;
    }

    // construct parameters from last track state
    return Result<void>::success();
  }

  /// Pointer to a logger that is owned by the parent, KalmanFilter
  std::shared_ptr<const Logger> m_logger{nullptr};

  /// Getter for the logger, to support logging macros
  const Logger& logger() const;
  
  private:

  template<unsigned int dimPrevTs, unsigned int dimTs>
  bool smooth(
  std::optional<std::error_code>& error,
  typename detail_lt::Types<dimTs, false>::CoefficientsMap tsFiltered,
  typename detail_lt::Types<dimTs, false>::CoefficientsMap tsSmoothed,
  typename detail_lt::Types<dimTs, false>::CovarianceMap tsFilteredCovariance,
  typename detail_lt::Types<dimTs, false>::CovarianceMap tsSmoothedCovariance,
  typename detail_lt::Types<dimPrevTs, false>::CoefficientsMap prevTsPredicted,  
  typename detail_lt::Types<dimPrevTs, false>::CoefficientsMap prevTsSmoothed,
  Eigen::Map<Eigen::Matrix<typename detail_lt::Types<dimPrevTs, false>::Scalar, dimPrevTs, dimTs, detail_lt::Types<dimPrevTs, false>::Flags>> prevTsJacobian,
  typename detail_lt::Types<dimPrevTs, false>::CovarianceMap prevTsPredictedCovariance,
  typename detail_lt::Types<dimPrevTs, false>::CovarianceMap prevTsSmoothedCovariance) const  {

		ACTS_VERBOSE("Calculate smoothing matrix:");
		ACTS_VERBOSE("Filtered covariance:\n" << tsFilteredCovariance);
		ACTS_VERBOSE("Jacobian:\n" << prevTsJacobian);
		ACTS_VERBOSE("Prev. predicted covariance\n"
					 << prevTsPredictedCovariance << "\n, inverse: \n"
					 << prevTsPredictedCovariance.inverse());
					 
		// Gain smoothing matrix
		// NB: The jacobian stored in a state is the jacobian from previous
		// state to this state in forward propagation
		const ActsMatrixD<dimTs, dimPrevTs> G = tsFilteredCovariance *
			prevTsJacobian.transpose() *
			prevTsPredictedCovariance.inverse();

		if (G.hasNaN()) {
		  error = KalmanFitterError::SmoothFailed;  // set to error
		  return false;                             // abort execution
		}

		ACTS_VERBOSE("Gain smoothing matrix G:\n" << G);					 
		ACTS_VERBOSE("Calculate smoothed parameters:");
		ACTS_VERBOSE("Filtered parameters: " << tsFiltered.transpose());
		ACTS_VERBOSE("Prev. smoothed parameters: "
					 << prevTsSmoothed.transpose());
		ACTS_VERBOSE("Prev. predicted parameters: "
					 << prevTsPredicted.transpose());
					 		
		// Calculate the smoothed parameters
		tsSmoothed = tsFiltered + G * (prevTsSmoothed - prevTsPredicted);

		ACTS_VERBOSE(
			"Smoothed parameters are: " << tsSmoothed.transpose());

		ACTS_VERBOSE("Calculate smoothed covariance:");
		ACTS_VERBOSE("Prev. smoothed covariance:\n"
					 << prevTsSmoothedCovariance);
					 
			// And the smoothed covariance
		tsSmoothedCovariance = tsFilteredCovariance - G * (prevTsPredictedCovariance -	prevTsSmoothedCovariance) * G.transpose();

		// Check if the covariance matrix is semi-positive definite.
		// If not, make one (could do more) attempt to replace it with the
		// nearest semi-positive def matrix,
		// but it could still be non semi-positive
		ActsSymMatrixD<dimTs> smoothedCov = tsSmoothedCovariance;
		if (not detail::covariance_helper<ActsSymMatrixD<dimTs>, 5>::validate(
				smoothedCov)) {
		  ACTS_DEBUG(
			  "Smoothed covariance is not positive definite. Could result in "
			  "negative covariance!");
		}
		// Reset smoothed covariance
		tsSmoothedCovariance = smoothedCov;
		ACTS_VERBOSE("Smoothed covariance is: \n"
					 << tsSmoothedCovariance);
		return true;
	}
};
}  // namespace Acts