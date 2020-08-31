// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

/// @brief void Measurement calibrator and converter
struct VoidKalmanComponents {
  /// @brief Public call mimicking a calibrator
  ///
  /// @tparam measurement_t Type of the measurement
  /// @tparam parameter_t Type of the parameters for calibration
  ///
  /// @param m Measurement to be moved through
  /// @param pars Parameters to be used for calibration
  ///
  /// @return void-calibrated measurement
  template <typename measurement_t, typename parameters_t>
  Result<measurement_t> operator()(measurement_t measurement,
                                   const parameters_t& /* parameters */) const {
    return measurement;
  }
};

/// @brief Void measurement calibrator for filtering
struct VoidMeasurementCalibrator {
  /// Main calibration call. In this implementation, it will dereference the
  /// given source link and expect it to result in something convertible to
  /// @c FittableMeasurement<source_link_t>.
  /// @tparam source_link_t Source link type which identifier the uncalibrated
  /// measurement
  /// @tparam parameters_t Parameters type (unused)
  /// @param sl Source link to turn into a measurement
  /// @param pars The parameters to calibrate with (unused)
  /// @note If the deref operator on @c source_link_t returns a reference, this
  /// will copy it before returning. If it is already returned by-value (for
  /// instance for a newly created measurement instance), return value
  /// optimizitaion should auto-move the result.
  /// @note This will not make the "calibrated" measurement point to the
  /// uncalibrated measurement via sourcelink, it's just a copy.
  template <typename source_link_t, typename parameters_t>
  FittableMeasurement<source_link_t> operator()(
      const source_link_t& sourceLink,
      const parameters_t& /* parameters */) const {
    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does fulfill SourceLinkConcept.");
    static_assert(
        concept ::converts_to<typename source_link_t::MeasurementType,
                              concept ::detail_slc::dereferenceable_t,
                              source_link_t>,
        "For DefaultMeasurementCalibrator, source link needs to implement "
        "dereference operator");
    return *sourceLink;
  }
};

/// @brief void Kalman updater
struct VoidKalmanUpdater {
  /// @brief Public call mimicking an updater
  ///
  /// @tparam measurement_t Type of the measurement to be used
  ///
  /// @param m The measurement
  /// @param predicted The predicted parameters
  ///
  /// @return The copied predicted parameters
  template <typename track_state_t>
  auto operator()(const GeometryContext& /*gctx*/, track_state_t& trackState,
                  const NavigationDirection& /*direction*/, bool boundState = true) const {
	if(boundState)
	{
		assert(trackState.hasBoundPredicted());
		assert(trackState.hasBoundFiltered());
      trackState.boundFiltered() = trackState.boundPredicted();
}
else
{
			assert(trackState.hasFreePredicted());
		assert(trackState.hasFreeFiltered());
	trackState.freeFiltered() = trackState.freePredicted();
}
  }
  /// Pointer to a logger that is owned by the parent, KalmanFilter
  std::shared_ptr<const Logger> m_logger{nullptr};
};

/// @brief void Kalman smoother
struct VoidKalmanSmoother {
  /// @brief Public call mimicking an smoother
  ///
  /// @tparam track_states_t Type of the track states
  ///
  /// @param states The track states to be smoothed
  ///
  /// @return The resulting
  template <typename parameters_t, typename track_states_t>
  const parameters_t* operator()(const GeometryContext& /*gctx*/, track_states_t& /* trackStates */, size_t /*entryIndex*/) const {
    return nullptr;
  }
  /// Pointer to a logger that is owned by the parent, KalmanFilter
  std::shared_ptr<const Logger> m_logger{nullptr};
};

/// @brief void outlier finder
struct VoidOutlierFinder {
  /// @brief Public call mimicking an outlier finder
  ///
  /// @tparam track_state_t Type of the track state
  ///
  /// @param trackState The trackState to investigate
  ///
  /// @return Whether it's outlier or not
  template <typename track_state_t>
  constexpr bool operator()(const track_state_t& /* trackState */) const {
    return false;
  }
};

}  // namespace Acts
