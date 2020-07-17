// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/StraightLineStepper.hpp"

#include "Acts/EventData/detail/coordinate_transformations.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"

namespace Acts {

auto StraightLineStepper::curvilinearState(State& state) const
    -> CurvilinearState {
  FreeVector parameters;
  parameters << state.pos[0], state.pos[1], state.pos[2], state.t, state.dir[0],
      state.dir[1], state.dir[2], (state.q != 0. ? state.q : 1.) / state.p;
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.jacDirToAngle, state.jacAngleToDir, parameters,
      state.covTransport, state.pathAccumulated);
}

auto StraightLineStepper::freeState(State& state) const -> FreeState {
  FreeVector parameters;
  parameters << state.pos[0], state.pos[1], state.pos[2], state.t, state.dir[0],
      state.dir[1], state.dir[2], (state.q != 0. ? state.q : 1.) / state.p;
  return detail::freeState(state.cov, state.jacobian, state.jacTransport,
                           state.derivative, state.jacToGlobal,
                           state.jacDirToAngle, state.jacAngleToDir, parameters,
                           state.covTransport, state.pathAccumulated);
}

auto StraightLineStepper::boundState(State& state, const Surface& surface) const
    -> BoundState {
  FreeVector parameters;
  parameters << state.pos[0], state.pos[1], state.pos[2], state.t, state.dir[0],
      state.dir[1], state.dir[2], (state.q != 0. ? state.q : 1.) / state.p;
  return detail::boundState(state.geoContext, state.cov, state.jacobian,
                            state.jacTransport, state.derivative,
                            state.jacToGlobal, state.jacDirToAngle,
                            state.jacAngleToDir, parameters, state.covTransport,
                            state.pathAccumulated, surface);
}

void StraightLineStepper::update(State& state, const FreeVector& parameters,
                                 const Covariance& covariance) const {
  state.pos = parameters.template segment<3>(eFreePos0);
  state.dir = parameters.template segment<3>(eFreeDir0).normalized();
  state.p = std::abs(1. / parameters[eFreeQOverP]);
  state.t = parameters[eFreeTime];

  state.cov = covariance;
}

void StraightLineStepper::update(State& state, const Vector3D& uposition,
                                 const Vector3D& udirection, double up,
                                 double time) const {
  state.pos = uposition;
  state.dir = udirection;
  state.p = up;
  state.t = time;
}

void StraightLineStepper::covarianceTransport(State& state,
                                              const Surface& surface) const {
  FreeVector parameters;
  parameters << state.pos[0], state.pos[1], state.pos[2], state.t, state.dir[0],
      state.dir[1], state.dir[2], (state.q != 0. ? state.q : 1.) / state.p;
  detail::covarianceTransport(state.geoContext, state.cov, state.jacobian,
                              state.jacTransport, state.derivative,
                              state.jacToGlobal, state.jacDirToAngle,
                              state.jacAngleToDir, parameters, surface);
}

void StraightLineStepper::resetState(State& state,
                                     const BoundVector& boundParams,
                                     const Covariance& cov,
                                     const Surface& surface,
                                     const NavigationDirection navDir,
                                     const double stepSize) const {
  using transformation = detail::coordinate_transformation;
  // Update the stepping state
  update(state,
         transformation::boundParameters2freeParameters(state.geoContext,
                                                        boundParams, surface),
         cov);
  state.navDir = navDir;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  state.jacToGlobal = BoundToFreeMatrix::Zero();
  surface.initJacobianToGlobal(state.geoContext, *state.jacToGlobal,
                               position(state), direction(state), boundParams);
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

void StraightLineStepper::resetState(State& state,
                                     const FreeVector& freeParams,
                                     const Covariance& cov,
                                     const NavigationDirection navDir,
                                     const double stepSize) const {
  // Update the stepping state
  update(state, freeParams, cov);
  state.navDir = navDir;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  state.jacToGlobal = std::nullopt;
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}
}  // namespace Acts
