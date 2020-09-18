// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Fitter/KalmanFitterError.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <variant>

namespace Acts {

/// @brief Update step of Kalman Filter using gain matrix formalism
class GainMatrixUpdater {
 public:
  /// Explicit constructor
  ///
  /// @param calibrator is the calibration struct/class that converts
  /// uncalibrated measurements into calibrated ones
  /// @param logger a logger instance
  GainMatrixUpdater(
      std::shared_ptr<const Logger> logger = std::shared_ptr<const Logger>(
          getDefaultLogger("GainMatrixUpdater", Logging::INFO).release()));

  /// @brief Public call operator for the boost visitor pattern
  ///
  /// @tparam track_state_t Type of the track state for the update
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param trackState the measured track state
  /// @param direction the navigation direction
  ///
  /// @return Bool indicating whether this update was 'successful'
  /// @note Non-'successful' updates could be holes or outliers,
  ///       which need to be treated differently in calling code.
  template <typename track_state_t>
  Result<void> operator()(
      const GeometryContext& /*gctx*/, track_state_t trackState,
      const NavigationDirection& direction = forward, bool boundState = true) const {
    ACTS_VERBOSE("Invoked GainMatrixUpdater");
    // let's make sure the types are consistent
    using SourceLink = typename track_state_t::SourceLink;
    using TrackStateProxy =
        typename MultiTrajectory<SourceLink>::TrackStateProxy;
    static_assert(std::is_same_v<track_state_t, TrackStateProxy>,
                  "Given track state type is not a track state proxy");

    // we should definitely have an uncalibrated measurement here
    assert(trackState.hasUncalibrated());
    // there should be a calibrated measurement
    assert(trackState.hasCalibrated());
    if(boundState)
    {
		// we should have predicted state set
		assert(trackState.hasBoundPredicted());
		// filtering should not have happened yet, but is allocated, therefore set
		assert(trackState.hasBoundFiltered());

		// read-only handles. Types are eigen maps to backing storage
		const auto predicted = trackState.boundPredicted();
		const auto predicted_covariance = trackState.boundPredictedCovariance();

		ACTS_VERBOSE("Predicted parameters: " << predicted.transpose());
		ACTS_VERBOSE("Predicted covariance:\n" << predicted_covariance);

		// read-write handles. Types are eigen maps into backing storage.
		// This writes directly into the trajectory storage
		auto filtered = trackState.boundFiltered();
		auto filtered_covariance = trackState.boundFilteredCovariance();

		return update<TrackStateProxy, 6>(trackState, predicted, predicted_covariance, filtered, filtered_covariance, direction);
	}
	else
	{
		// we should have predicted state set
		assert(trackState.hasFreePredicted());
		// filtering should not have happened yet, but is allocated, therefore set
		assert(trackState.hasFreeFiltered());

		// read-only handles. Types are eigen maps to backing storage
		const auto predicted = trackState.freePredicted();
		const auto predicted_covariance = trackState.freePredictedCovariance();

		ACTS_VERBOSE("Predicted parameters: " << predicted.transpose());
		ACTS_VERBOSE("Predicted covariance:\n" << predicted_covariance);

		// read-write handles. Types are eigen maps into backing storage.
		// This writes directly into the trajectory storage
		auto filtered = trackState.freeFiltered();
		auto filtered_covariance = trackState.freeFilteredCovariance();

		return update<TrackStateProxy, 8>(trackState, predicted, predicted_covariance, filtered, filtered_covariance, direction);
  }
}

  /// Pointer to a logger that is owned by the parent, KalmanFilter
  std::shared_ptr<const Logger> m_logger{nullptr};

  /// Getter for the logger, to support logging macros
  const Logger& logger() const;
  
private:
  /// @brief This function runs the actual gain matrix update
  ///
  /// @tparam track_state_t Type of the track state
  /// @tparam parameter_size_t Size of the parameter vectors
  ///
  /// @param [in] trackState Current track state
  /// @param [in] predicted Predicted parameter vector
  /// @param [in] predictedCovariance Predicted covariance matrix
  /// @param [in, out] filtered Filtered parameter vector
  /// @param [in, out] filteredCovariance Filtered covariance matrix
  /// @param [in] direction Navigation direction
  ///
  /// @return Bool indicating whether this update was 'successful'
  template <typename track_state_t, int parameter_size_t>
  Result<void> update(
      track_state_t& trackState, 
      typename detail_lt::Types<parameter_size_t, false>::CoefficientsMap predicted, 
      typename detail_lt::Types<parameter_size_t, false>::CovarianceMap predictedCovariance, 
      typename detail_lt::Types<parameter_size_t, false>::CoefficientsMap filtered, 
      typename detail_lt::Types<parameter_size_t, false>::CovarianceMap filteredCovariance,
      const NavigationDirection& direction = forward ) const {
	std::optional<std::error_code> error{std::nullopt};  // assume ok
	visit_measurement(
			trackState.calibrated(), trackState.calibratedCovariance(),
			trackState.calibratedSize(),
			[&](const auto calibrated, const auto calibrated_covariance) {
			  constexpr size_t measdim = decltype(calibrated)::RowsAtCompileTime;
			  using cov_t = ActsSymMatrixD<measdim>;
			  using par_t = ActsVectorD<measdim>;

			  ACTS_VERBOSE("Measurement dimension: " << measdim);
			  ACTS_VERBOSE("Calibrated measurement: " << calibrated.transpose());
			  ACTS_VERBOSE("Calibrated measurement covariance:\n"
						   << calibrated_covariance);

			  const ActsMatrixD<measdim, parameter_size_t> H =
				  trackState.projector()
					  .template topLeftCorner<measdim, parameter_size_t>();

			  ACTS_VERBOSE("Measurement projector H:\n" << H);
			  ActsMatrixD<parameter_size_t, measdim> K;
			  ActsMatrixD<parameter_size_t, measdim> Q;
			  ActsMatrixD<parameter_size_t, parameter_size_t> B;
			  if constexpr (parameter_size_t == 8 && measdim == 3)
			  {
				  Vector3D dir = predicted.template segment<3>(eFreeDir0);
				  ActsSymMatrixD<3> dirMat = dir * dir.transpose();
				  ActsSymMatrixD<3> dirMat2 = dirMat + predictedCovariance.template block<3, 3>(eFreeDir0, eFreeDir0);
				  double v = calibrated_covariance.diagonal().sum();
				  B = H.transpose() * v * dirMat2 * H;
				  predictedCovariance = predictedCovariance + B;
				  
				  //~ K =
				  //~ predictedCovariance * H.transpose() *
				  //~ (H * predictedCovariance * H.transpose() + calibrated_covariance)
					  //~ .inverse();
				  //~ K =
				  //~ (predictedCovariance * H.transpose() - Q) *
				  //~ (H * predictedCovariance * H.transpose() + calibrated_covariance - H * Q - Q.transpose() * H.transpose())
					  //~ .inverse();
				  //~ K =
				  //~ (predictedCovariance - B) * H.transpose() *
				  //~ (H * predictedCovariance * H.transpose() + calibrated_covariance - H * Q - Q.transpose() * H.transpose())
					  //~ .inverse();
				  K =
				  (predictedCovariance - B) * H.transpose() *
				  (H * predictedCovariance * H.transpose() + calibrated_covariance)
					  .inverse();
			  } else {
				K =
				  predictedCovariance * H.transpose() *
				  (H * predictedCovariance * H.transpose() + calibrated_covariance)
					  .inverse();
			}

			  ACTS_VERBOSE("Gain Matrix K:\n" << K);

			  if (K.hasNaN()) {
				error =
					(direction == forward)
						? KalmanFitterError::ForwardUpdateFailed
						: KalmanFitterError::BackwardUpdateFailed;  // set to error
				return false;  // abort execution
			  }
			  filtered = predicted + K * (calibrated - H * predicted);
			  if constexpr (parameter_size_t == 8 && measdim == 3)
			  {
				  //~ filteredCovariance = (ActsSymMatrixD<parameter_size_t>::Identity() - K * H) * predictedCovariance + K * H * B;
				  //~ filteredCovariance = (ActsSymMatrixD<parameter_size_t>::Identity() - K * H) * predictedCovariance + 2. * B - B * H.transpose() * K.transpose();
				  filteredCovariance = (ActsSymMatrixD<parameter_size_t>::Identity() - K * H) * (predictedCovariance - B);
				  predictedCovariance = predictedCovariance - B;
//~ std::cout << "(ActsSymMatrixD<parameter_size_t>::Identity() - K * H)\n" << (ActsSymMatrixD<parameter_size_t>::Identity() - K * H) << std::endl;
//~ std::cout << "Predicted:\n" << predictedCovariance << std::endl;
//~ std::cout << "K * Q.transpose()\n" << K * Q.transpose() << std::endl;
//~ std::cout << "Filtered:\n" << filteredCovariance << std::endl;
			  } else {
				  filteredCovariance =
					  (ActsSymMatrixD<parameter_size_t>::Identity() - K * H) * predictedCovariance;
					 }
			  ACTS_VERBOSE("Filtered parameters: " << filtered.transpose());
			  ACTS_VERBOSE("Filtered covariance:\n" << filteredCovariance);

			  // calculate filtered residual
			  par_t residual(trackState.calibratedSize());
			  residual = (calibrated - H * filtered);
			  ACTS_VERBOSE("Residual: " << residual.transpose());

			  trackState.chi2() =
				  (residual.transpose() *
				   ((cov_t::Identity() - H * K) * calibrated_covariance).inverse() *
				   residual)
					  .value();
// TODO: If a direction component was changed, a normalisation is required
			  ACTS_VERBOSE("Chi2: " << trackState.chi2());
			  return true;  // continue execution
			});
	if (error) {
      // error is set, return result
      return *error;
    }
    // always succeed, no outlier logic yet
    return Result<void>::success();
  }
};

}  // namespace Acts
