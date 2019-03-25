// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Units.hpp"

// This is based original stepper code from the ATLAS RungeKuttePropagagor
namespace Acts {

/// @brief the AtlasSTEPStepper implementation for the
template <typename bfield_t>
class AtlasSTEPStepper
{

public:
  using Jacobian = ActsMatrixD<5, 5>;
  using Cstep    = detail::ConstrainedStep;

  /// @brief Nested State struct for the local caching
  struct State
  {

    /// Constructor
    ///
    /// @tparams Type of TrackParameters
    ///
    /// @param[in] pars Input parameters
    /// @param[in] ndir The navigation direction w.r.t. parameters
    /// @param[in] ssize the steps size limitation
    template <typename Parameters>
    State(const Parameters&   pars,
          NavigationDirection ndir  = forward,
          double              ssize = std::numeric_limits<double>::max())
      : state_ready(false)
      , navDir(ndir)
      , useJacobian(false)
      , step(0.)
      , maxPathLength(0.)
      , mcondition(false)
      , needgradient(false)
      , newfield(true)
      , field(0., 0., 0.)
      , covariance(nullptr)
      , stepSize(ndir * std::abs(ssize))
    {
      // The rest of this constructor is copy&paste of AtlasSTEPStepper::update() -
      // this is a nasty but working solution for the stepper state without
      // functions

      const ActsVectorD<3> pos = pars.position();
      const auto           Vp  = pars.parameters();

      double Sf, Cf, Ce, Se;
      Sf = sin(Vp(2));
      Cf = cos(Vp(2));
      Se = sin(Vp(3));
      Ce = cos(Vp(3));

      pVector[0] = pos(0);
      pVector[1] = pos(1);
      pVector[2] = pos(2);
      pVector[3] = Cf * Se;
      pVector[4] = Sf * Se;
      pVector[5] = Ce;
      pVector[6] = Vp[4];

      // @todo: remove magic numbers - is that the charge ?
      if (std::abs(pVector[6]) < .000000000000001) {
        pVector[6] < 0. ? pVector[6] = -.000000000000001
                        : pVector[6] = .000000000000001;
      }

      // prepare the jacobian if we have a covariance
      if (pars.covariance()) {
        // copy the covariance matrix
        covariance  = new ActsSymMatrixD<NGlobalPars>(*pars.covariance());
        useJacobian = true;
        const auto transform = pars.referenceFrame();

        pVector[7]  = transform(0, eLOC_0);
        pVector[14] = transform(0, eLOC_1);
        pVector[21] = 0.;
        pVector[28] = 0.;
        pVector[35] = 0.;  // dX /

        pVector[8]  = transform(1, eLOC_0);
        pVector[15] = transform(1, eLOC_1);
        pVector[22] = 0.;
        pVector[29] = 0.;
        pVector[36] = 0.;  // dY /

        pVector[9]  = transform(2, eLOC_0);
        pVector[16] = transform(2, eLOC_1);
        pVector[23] = 0.;
        pVector[30] = 0.;
        pVector[37] = 0.;  // dZ /

        pVector[10] = 0.;
        pVector[17] = 0.;
        pVector[24] = -Sf * Se;  // - sin(phi) * cos(theta)
        pVector[31] = Cf * Ce;   // cos(phi) * cos(theta)
        pVector[38] = 0.;        // dAx/

        pVector[11] = 0.;
        pVector[18] = 0.;
        pVector[25] = Cf * Se;  // cos(phi) * sin(theta)
        pVector[32] = Sf * Ce;  // sin(phi) * cos(theta)
        pVector[39] = 0.;       // dAy/

        pVector[12] = 0.;
        pVector[19] = 0.;
        pVector[26] = 0.;
        pVector[33] = -Se;  // - sin(theta)
        pVector[40] = 0.;   // dAz/

        pVector[13] = 0.;
        pVector[20] = 0.;
        pVector[27] = 0.;
        pVector[34] = 0.;
        pVector[41] = 1.;  // dCM/

        pVector[42] = 0.;
        pVector[43] = 0.;
        pVector[44] = 0.;

        // special treatment for surface types
        const auto& surface = pars.referenceSurface();
        // the disc needs polar coordinate adaptations
        if (surface.type() == Surface::Disc) {
          double lCf   = cos(Vp[1]);
          double lSf   = sin(Vp[1]);
          double Ax[3] = {transform(0, 0), transform(1, 0), transform(2, 0)};
          double Ay[3] = {transform(0, 1), transform(1, 1), transform(2, 1)};
          double d0    = lCf * Ax[0] + lSf * Ay[0];
          double d1    = lCf * Ax[1] + lSf * Ay[1];
          double d2    = lCf * Ax[2] + lSf * Ay[2];
          pVector[7]   = d0;
          pVector[8]   = d1;
          pVector[9]   = d2;
          pVector[14]  = Vp[0] * (lCf * Ay[0] - lSf * Ax[0]);
          pVector[15]  = Vp[0] * (lCf * Ay[1] - lSf * Ax[1]);
          pVector[16]  = Vp[0] * (lCf * Ay[2] - lSf * Ax[2]);
        }
        // the line needs components that relate direction change
        // with global frame change
        if (surface.type() == Surface::Perigee
            || surface.type() == Surface::Straw) {

          // sticking to the nomenclature of the original RkPropagator
          // - axis pointing along the drift/transverse direction
          double B[3] = {transform(0, 0), transform(1, 0), transform(2, 0)};
          // - axis along the straw
          double A[3] = {transform(0, 1), transform(1, 1), transform(2, 1)};
          // - normal vector of the reference frame
          double C[3] = {transform(0, 2), transform(1, 2), transform(2, 2)};

          // projection of direction onto normal vector of reference frame
          double PC = pVector[3] * C[0] + pVector[4] * C[1] + pVector[5] * C[2];
          double Bn = 1. / PC;

          double Bx2 = -A[2] * pVector[25];
          double Bx3 = A[1] * pVector[33] - A[2] * pVector[32];

          double By2 = A[2] * pVector[24];
          double By3 = A[2] * pVector[31] - A[0] * pVector[33];

          double Bz2 = A[0] * pVector[25] - A[1] * pVector[24];
          double Bz3 = A[0] * pVector[32] - A[1] * pVector[31];

          double B2 = B[0] * Bx2 + B[1] * By2 + B[2] * Bz2;
          double B3 = B[0] * Bx3 + B[1] * By3 + B[2] * Bz3;

          Bx2 = (Bx2 - B[0] * B2) * Bn;
          Bx3 = (Bx3 - B[0] * B3) * Bn;
          By2 = (By2 - B[1] * B2) * Bn;
          By3 = (By3 - B[1] * B3) * Bn;
          Bz2 = (Bz2 - B[2] * B2) * Bn;
          Bz3 = (Bz3 - B[2] * B3) * Bn;

          //  /dPhi      |     /dThe       |
          pVector[21] = Bx2 * Vp[0];
          pVector[28] = Bx3 * Vp[0];  // dX/
          pVector[22] = By2 * Vp[0];
          pVector[29] = By3 * Vp[0];  // dY/
          pVector[23] = Bz2 * Vp[0];
          pVector[30] = Bz3 * Vp[0];  // dZ/
        }
      }
      // now declare the state as ready
      state_ready = true;
    }

    // optimisation that init is not called twice
    bool state_ready = false;
    // configuration
    NavigationDirection navDir;
    bool                useJacobian;
    double              step;
    double              maxPathLength;
    bool                mcondition;
    bool                needgradient;
    bool                newfield;
    // internal parameters to be used
    Vector3D field;
    double   pVector[64];
    // result
    double parameters[NGlobalPars] = {0., 0., 0., 0., 0.};
    const ActsSymMatrixD<NGlobalPars>* covariance;
    double                             jacobian[NGlobalPars * NGlobalPars];
    
    bool m_energyLoss = true;
    bool errorPropagation = true;
    bool firstStep = true;
    ActsVectorF<12> BG1;
    double m_momentumCutOff = 0.;
    bool m_includeGgradient = false;
    bool m_includeBgradients = false;
    double m_delIoni, m_delRad, m_kazL;
    bool m_MPV = false;
    double m_sigmaIoni;
    double m_sigmaRad;

    /// Lazily initialized cache for the magnetic field
    /// It caches the current magnetic field cell and stays (and interpolates)
    ///  within as long as this is valid. See step() code for details.
    typename bfield_t::Cache fieldCache{};

    // accummulated path length cache
    double pathAccumulated = 0.;

    // adaptive step size of the runge-kutta integration
    Cstep stepSize = std::numeric_limits<double>::max();

    /// Debug output
    /// the string where debug messages are stored (optionally)
    bool        debug       = false;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    size_t debugPfxWidth = 30;
    size_t debugMsgWidth = 50;
  };

  Vector3D
  position(const State& state) const
  {
    return Vector3D(state.pVector[0], state.pVector[1], state.pVector[2]);
  }

  Vector3D
  direction(const State& state) const
  {
    return Vector3D(state.pVector[3], state.pVector[4], state.pVector[5]);
  }

  double
  momentum(const State& state) const
  {
    return 1. / std::abs(state.pVector[6]);
  }

  /// Charge access
  double
  charge(const State& state) const
  {
    return state.pVector[6] > 0. ? 1. : -1.;
  }

  /// Method to update momentum, direction and p
  ///
  /// @param uposition the updated position
  /// @param udirection the updated direction
  /// @param p the updated momentum value
  void
  update(State&          state,
         const Vector3D& uposition,
         const Vector3D& udirection,
         double          up)
  {
    // update the vector
    state.pVector[0] = uposition[0];
    state.pVector[1] = uposition[1];
    state.pVector[2] = uposition[2];
    state.pVector[3] = udirection[0];
    state.pVector[4] = udirection[1];
    state.pVector[5] = udirection[2];
    state.pVector[6] = charge(state) / up;
  }

  /// Return a corrector
  VoidIntersectionCorrector
  corrector(State& /*unused*/) const
  {
    return VoidIntersectionCorrector();
  }

  /// The state update method
  ///
  /// @param [in] pars The new track parameters at start
  template <typename Parameters>
  void
  update(State& state, const Parameters& pars)
  {
    // state is ready - noting to do
    if (state.state_ready) {
      return;
    }

    const ActsVectorD<3> pos = pars.position();
    const auto           Vp  = pars.parameters();

    double Sf, Cf, Ce, Se;
    Sf = sin(Vp(2));
    Cf = cos(Vp(2));
    Se = sin(Vp(3));
    Ce = cos(Vp(3));

    state.pVector[0] = pos(0);
    state.pVector[1] = pos(1);
    state.pVector[2] = pos(2);
    state.pVector[3] = Cf * Se;
    state.pVector[4] = Sf * Se;
    state.pVector[5] = Ce;
    state.pVector[6] = Vp[4];

    // @todo: remove magic numbers - is that the charge ?
    if (std::abs(state.pVector[6]) < .000000000000001) {
      state.pVector[6] < 0. ? state.pVector[6] = -.000000000000001
                            : state.pVector[6] = .000000000000001;
    }

    // prepare the jacobian if we have a covariance
    if (pars.covariance()) {
      // copy the covariance matrix
      state.covariance  = new ActsSymMatrixD<NGlobalPars>(*pars.covariance());
      state.useJacobian = true;
      const auto transform = pars.referenceFrame();

      state.pVector[7]  = transform(0, eLOC_0);
      state.pVector[14] = transform(0, eLOC_1);
      state.pVector[21] = 0.;
      state.pVector[28] = 0.;
      state.pVector[35] = 0.;  // dX /

      state.pVector[8]  = transform(1, eLOC_0);
      state.pVector[15] = transform(1, eLOC_1);
      state.pVector[22] = 0.;
      state.pVector[29] = 0.;
      state.pVector[36] = 0.;  // dY /

      state.pVector[9]  = transform(2, eLOC_0);
      state.pVector[16] = transform(2, eLOC_1);
      state.pVector[23] = 0.;
      state.pVector[30] = 0.;
      state.pVector[37] = 0.;  // dZ /

      state.pVector[10] = 0.;
      state.pVector[17] = 0.;
      state.pVector[24] = -Sf * Se;  // - sin(phi) * cos(theta)
      state.pVector[31] = Cf * Ce;   // cos(phi) * cos(theta)
      state.pVector[38] = 0.;        // dAx/

      state.pVector[11] = 0.;
      state.pVector[18] = 0.;
      state.pVector[25] = Cf * Se;  // cos(phi) * sin(theta)
      state.pVector[32] = Sf * Ce;  // sin(phi) * cos(theta)
      state.pVector[39] = 0.;       // dAy/

      state.pVector[12] = 0.;
      state.pVector[19] = 0.;
      state.pVector[26] = 0.;
      state.pVector[33] = -Se;  // - sin(theta)
      state.pVector[40] = 0.;   // dAz/

      state.pVector[13] = 0.;
      state.pVector[20] = 0.;
      state.pVector[27] = 0.;
      state.pVector[34] = 0.;
      state.pVector[41] = 1.;  // dCM/

      state.pVector[42] = 0.;
      state.pVector[43] = 0.;
      state.pVector[44] = 0.;

      // special treatment for surface types
      const auto& surface = pars.referenceSurface();
      // the disc needs polar coordinate adaptations
      if (surface.type() == Surface::Disc) {
        double lCf        = cos(Vp[1]);
        double lSf        = sin(Vp[1]);
        double Ax[3]      = {transform(0, 0), transform(1, 0), transform(2, 0)};
        double Ay[3]      = {transform(0, 1), transform(1, 1), transform(2, 1)};
        double d0         = lCf * Ax[0] + lSf * Ay[0];
        double d1         = lCf * Ax[1] + lSf * Ay[1];
        double d2         = lCf * Ax[2] + lSf * Ay[2];
        state.pVector[7]  = d0;
        state.pVector[8]  = d1;
        state.pVector[9]  = d2;
        state.pVector[14] = Vp[0] * (lCf * Ay[0] - lSf * Ax[0]);
        state.pVector[15] = Vp[0] * (lCf * Ay[1] - lSf * Ax[1]);
        state.pVector[16] = Vp[0] * (lCf * Ay[2] - lSf * Ax[2]);
      }
      // the line needs components that relate direction change
      // with global frame change
      if (surface.type() == Surface::Perigee
          || surface.type() == Surface::Straw) {

        // sticking to the nomenclature of the original RkPropagator
        // - axis pointing along the drift/transverse direction
        double B[3] = {transform(0, 0), transform(1, 0), transform(2, 0)};
        // - axis along the straw
        double A[3] = {transform(0, 1), transform(1, 1), transform(2, 1)};
        // - normal vector of the reference frame
        double C[3] = {transform(0, 2), transform(1, 2), transform(2, 2)};

        // projection of direction onto normal vector of reference frame
        double PC = state.pVector[3] * C[0] + state.pVector[4] * C[1]
            + state.pVector[5] * C[2];
        double Bn = 1. / PC;

        double Bx2 = -A[2] * state.pVector[25];
        double Bx3 = A[1] * state.pVector[33] - A[2] * state.pVector[32];

        double By2 = A[2] * state.pVector[24];
        double By3 = A[2] * state.pVector[31] - A[0] * state.pVector[33];

        double Bz2 = A[0] * state.pVector[25] - A[1] * state.pVector[24];
        double Bz3 = A[0] * state.pVector[32] - A[1] * state.pVector[31];

        double B2 = B[0] * Bx2 + B[1] * By2 + B[2] * Bz2;
        double B3 = B[0] * Bx3 + B[1] * By3 + B[2] * Bz3;

        Bx2 = (Bx2 - B[0] * B2) * Bn;
        Bx3 = (Bx3 - B[0] * B3) * Bn;
        By2 = (By2 - B[1] * B2) * Bn;
        By3 = (By3 - B[1] * B3) * Bn;
        Bz2 = (Bz2 - B[2] * B2) * Bn;
        Bz3 = (Bz3 - B[2] * B3) * Bn;

        //  /dPhi      |     /dThe       |
        state.pVector[21] = Bx2 * Vp[0];
        state.pVector[28] = Bx3 * Vp[0];  // dX/
        state.pVector[22] = By2 * Vp[0];
        state.pVector[29] = By3 * Vp[0];  // dY/
        state.pVector[23] = Bz2 * Vp[0];
        state.pVector[30] = Bz3 * Vp[0];  // dZ/
      }
    }
    // now declare the state as ready
    state.state_ready = true;
  }

  template <typename T, typename S = int>
  using state_type = State;

  template <typename T>
  using step_parameter_type = CurvilinearParameters;

  // This struct is a meta-function which normally maps to BoundParameters...
  template <typename T, typename S>
  struct s
  {
    using type = BoundParameters;
  };

  // Unless S is int, then it maps to CurvilinearParameters ...
  template <typename T>
  struct s<T, int>
  {
    using type = CurvilinearParameters;
  };

  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;

  /// Convert the propagation state (global) to curvilinear parameters
  /// This is called by the propagator
  ///
  /// @tparam result_t Type of the propagator result to be filled
  ///
  /// @param[in,out] state The stepper state
  /// @param[in,out] result The propagator result object to be filled
  template <typename result_t>
  void
  convert(State& state, result_t& result) const
  {
    // the convert method invalidates the state (in case it's reused)
    state.state_ready = false;
    //
    Acts::Vector3D gp(state.pVector[0], state.pVector[1], state.pVector[2]);
    Acts::Vector3D mom(state.pVector[3], state.pVector[4], state.pVector[5]);
    mom /= std::abs(state.pVector[6]);

    double P[45];
    for (unsigned int i = 0; i < 45; ++i) {
      P[i] = state.pVector[i];
    }

    std::unique_ptr<const ActsSymMatrixD<NGlobalPars>> cov = nullptr;
    if (state.covariance) {
      double p = 1. / P[6];
      P[35] *= p;
      P[36] *= p;
      P[37] *= p;
      P[38] *= p;
      P[39] *= p;
      P[40] *= p;

      double An = sqrt(P[3] * P[3] + P[4] * P[4]);
      double Ax[3];
      if (An != 0.) {
        Ax[0] = -P[4] / An;
        Ax[1] = P[3] / An;
        Ax[2] = 0.;
      } else {
        Ax[0] = 1.;
        Ax[1] = 0.;
        Ax[2] = 0.;
      }

      double Ay[3] = {-Ax[1] * P[5], Ax[0] * P[5], An};
      double S[3]  = {P[3], P[4], P[5]};

      double A = P[3] * S[0] + P[4] * S[1] + P[5] * S[2];
      if (A != 0.) {
        A = 1. / A;
      }
      S[0] *= A;
      S[1] *= A;
      S[2] *= A;

      double s0 = P[7] * S[0] + P[8] * S[1] + P[9] * S[2];
      double s1 = P[14] * S[0] + P[15] * S[1] + P[16] * S[2];
      double s2 = P[21] * S[0] + P[22] * S[1] + P[23] * S[2];
      double s3 = P[28] * S[0] + P[29] * S[1] + P[30] * S[2];
      double s4 = P[35] * S[0] + P[36] * S[1] + P[37] * S[2];

      P[7] -= (s0 * P[3]);
      P[8] -= (s0 * P[4]);
      P[9] -= (s0 * P[5]);
      P[10] -= (s0 * P[42]);
      P[11] -= (s0 * P[43]);
      P[12] -= (s0 * P[44]);
      P[14] -= (s1 * P[3]);
      P[15] -= (s1 * P[4]);
      P[16] -= (s1 * P[5]);
      P[17] -= (s1 * P[42]);
      P[18] -= (s1 * P[43]);
      P[19] -= (s1 * P[44]);
      P[21] -= (s2 * P[3]);
      P[22] -= (s2 * P[4]);
      P[23] -= (s2 * P[5]);
      P[24] -= (s2 * P[42]);
      P[25] -= (s2 * P[43]);
      P[26] -= (s2 * P[44]);
      P[28] -= (s3 * P[3]);
      P[29] -= (s3 * P[4]);
      P[30] -= (s3 * P[5]);
      P[31] -= (s3 * P[42]);
      P[32] -= (s3 * P[43]);
      P[33] -= (s3 * P[44]);
      P[35] -= (s4 * P[3]);
      P[36] -= (s4 * P[4]);
      P[37] -= (s4 * P[5]);
      P[38] -= (s4 * P[42]);
      P[39] -= (s4 * P[43]);
      P[40] -= (s4 * P[44]);

      double P3, P4, C = P[3] * P[3] + P[4] * P[4];
      if (C > 1.e-20) {
        C  = 1. / C;
        P3 = P[3] * C;
        P4 = P[4] * C;
        C  = -sqrt(C);
      } else {
        C  = -1.e10;
        P3 = 1.;
        P4 = 0.;
      }

      // Jacobian production
      //
      state.jacobian[0] = Ax[0] * P[7] + Ax[1] * P[8];    // dL0/dL0
      state.jacobian[1] = Ax[0] * P[14] + Ax[1] * P[15];  // dL0/dL1
      state.jacobian[2] = Ax[0] * P[21] + Ax[1] * P[22];  // dL0/dPhi
      state.jacobian[3] = Ax[0] * P[28] + Ax[1] * P[29];  // dL0/dThe
      state.jacobian[4] = Ax[0] * P[35] + Ax[1] * P[36];  // dL0/dCM

      state.jacobian[5]
          = Ay[0] * P[7] + Ay[1] * P[8] + Ay[2] * P[9];  // dL1/dL0
      state.jacobian[6]
          = Ay[0] * P[14] + Ay[1] * P[15] + Ay[2] * P[16];  // dL1/dL1
      state.jacobian[7]
          = Ay[0] * P[21] + Ay[1] * P[22] + Ay[2] * P[23];  // dL1/dPhi
      state.jacobian[8]
          = Ay[0] * P[28] + Ay[1] * P[29] + Ay[2] * P[30];  // dL1/dThe
      state.jacobian[9]
          = Ay[0] * P[35] + Ay[1] * P[36] + Ay[2] * P[37];  // dL1/dCM

      state.jacobian[10] = P3 * P[11] - P4 * P[10];  // dPhi/dL0
      state.jacobian[11] = P3 * P[18] - P4 * P[17];  // dPhi/dL1
      state.jacobian[12] = P3 * P[25] - P4 * P[24];  // dPhi/dPhi
      state.jacobian[13] = P3 * P[32] - P4 * P[31];  // dPhi/dThe
      state.jacobian[14] = P3 * P[39] - P4 * P[38];  // dPhi/dCM

      state.jacobian[15] = C * P[12];  // dThe/dL0
      state.jacobian[16] = C * P[19];  // dThe/dL1
      state.jacobian[17] = C * P[26];  // dThe/dPhi
      state.jacobian[18] = C * P[33];  // dThe/dThe
      state.jacobian[19] = C * P[40];  // dThe/dCM

      state.jacobian[20] = 0.;     // dCM /dL0
      state.jacobian[21] = 0.;     // dCM /dL1
      state.jacobian[22] = 0.;     // dCM /dPhi
      state.jacobian[23] = 0.;     // dCM /dTheta
      state.jacobian[24] = P[41];  // dCM /dCM

      Eigen::
          Map<Eigen::Matrix<double, NGlobalPars, NGlobalPars, Eigen::RowMajor>>
              J(state.jacobian);

      cov = std::make_unique<const ActsSymMatrixD<NGlobalPars>>(
          J * (*state.covariance) * J.transpose());
      // Optionally : fill the jacobian
      result.transportJacobian = std::make_unique<const Jacobian>(std::move(J));
    }

    // Fill the result
    result.endParameters = std::make_unique<const CurvilinearParameters>(
        std::move(cov), gp, mom, charge(state));
  }

  /// Convert the propagation state to track parameters at a certain surface
  ///
  /// @tparam result_t Type of the propagator result to be filled
  /// @tparam surface_t Type of the surface
  ///
  /// @param [in,out] state Propagation state used
  /// @param [in,out] result Result object from the propagator
  /// @param [in] s Destination surface to which the conversion is done
  template <typename result_t, typename surface_t>
  void
  convert(State& state, result_t& result, const surface_t& surface) const
  {

    // the convert method invalidates the state (in case it's reused)
    state.state_ready = false;

    /// The transport of the position
    Acts::Vector3D gp(state.pVector[0], state.pVector[1], state.pVector[2]);
    Acts::Vector3D mom(state.pVector[3], state.pVector[4], state.pVector[5]);
    mom /= std::abs(state.pVector[6]);

    // The transport of the covariance
    std::unique_ptr<const ActsSymMatrixD<5>> cov = nullptr;
    if (state.covariance) {
      double p = 1. / state.pVector[6];
      state.pVector[35] *= p;
      state.pVector[36] *= p;
      state.pVector[37] *= p;
      state.pVector[38] *= p;
      state.pVector[39] *= p;
      state.pVector[40] *= p;

      const auto fFrame = surface.referenceFrame(gp, mom);

      double Ax[3] = {fFrame(0, 0), fFrame(1, 0), fFrame(2, 0)};
      double Ay[3] = {fFrame(0, 1), fFrame(1, 1), fFrame(2, 1)};
      double S[3]  = {fFrame(0, 2), fFrame(1, 2), fFrame(2, 2)};

      // this is the projection of direction onto the local normal vector
      double A = state.pVector[3] * S[0] + state.pVector[4] * S[1]
          + state.pVector[5] * S[2];

      if (A != 0.) {
        A = 1. / A;
      }

      S[0] *= A;
      S[1] *= A;
      S[2] *= A;

      double s0 = state.pVector[7] * S[0] + state.pVector[8] * S[1]
          + state.pVector[9] * S[2];
      double s1 = state.pVector[14] * S[0] + state.pVector[15] * S[1]
          + state.pVector[16] * S[2];
      double s2 = state.pVector[21] * S[0] + state.pVector[22] * S[1]
          + state.pVector[23] * S[2];
      double s3 = state.pVector[28] * S[0] + state.pVector[29] * S[1]
          + state.pVector[30] * S[2];
      double s4 = state.pVector[35] * S[0] + state.pVector[36] * S[1]
          + state.pVector[37] * S[2];

      // in case of line-type surfaces - we need to take into account that
      // the reference frame changes with variations of all local
      // parameters
      if (surface.type() == Surface::Straw
          || surface.type() == Surface::Perigee) {
        // vector from position to center
        double x = state.pVector[0] - surface.center().x();
        double y = state.pVector[1] - surface.center().y();
        double z = state.pVector[2] - surface.center().z();

        // this is the projection of the direction onto the local y axis
        double d = state.pVector[3] * Ay[0] + state.pVector[4] * Ay[1]
            + state.pVector[5] * Ay[2];

        // this is cos(beta)
        double a = (1. - d) * (1. + d);
        if (a != 0.) {
          a = 1. / a;  // i.e. 1./(1-d^2)
        }

        // that's the modified norm vector
        double X = d * Ay[0] - state.pVector[3];  //
        double Y = d * Ay[1] - state.pVector[4];  //
        double Z = d * Ay[2] - state.pVector[5];  //

        // d0 to d1
        double d0 = state.pVector[10] * Ay[0] + state.pVector[11] * Ay[1]
            + state.pVector[12] * Ay[2];
        double d1 = state.pVector[17] * Ay[0] + state.pVector[18] * Ay[1]
            + state.pVector[19] * Ay[2];
        double d2 = state.pVector[24] * Ay[0] + state.pVector[25] * Ay[1]
            + state.pVector[26] * Ay[2];
        double d3 = state.pVector[31] * Ay[0] + state.pVector[32] * Ay[1]
            + state.pVector[33] * Ay[2];
        double d4 = state.pVector[38] * Ay[0] + state.pVector[39] * Ay[1]
            + state.pVector[40] * Ay[2];

        s0 = (((state.pVector[7] * X + state.pVector[8] * Y
                + state.pVector[9] * Z)
               + x * (d0 * Ay[0] - state.pVector[10]))
              + (y * (d0 * Ay[1] - state.pVector[11])
                 + z * (d0 * Ay[2] - state.pVector[12])))
            * (-a);

        s1 = (((state.pVector[14] * X + state.pVector[15] * Y
                + state.pVector[16] * Z)
               + x * (d1 * Ay[0] - state.pVector[17]))
              + (y * (d1 * Ay[1] - state.pVector[18])
                 + z * (d1 * Ay[2] - state.pVector[19])))
            * (-a);
        s2 = (((state.pVector[21] * X + state.pVector[22] * Y
                + state.pVector[23] * Z)
               + x * (d2 * Ay[0] - state.pVector[24]))
              + (y * (d2 * Ay[1] - state.pVector[25])
                 + z * (d2 * Ay[2] - state.pVector[26])))
            * (-a);
        s3 = (((state.pVector[28] * X + state.pVector[29] * Y
                + state.pVector[30] * Z)
               + x * (d3 * Ay[0] - state.pVector[31]))
              + (y * (d3 * Ay[1] - state.pVector[32])
                 + z * (d3 * Ay[2] - state.pVector[33])))
            * (-a);
        s4 = (((state.pVector[35] * X + state.pVector[36] * Y
                + state.pVector[37] * Z)
               + x * (d4 * Ay[0] - state.pVector[38]))
              + (y * (d4 * Ay[1] - state.pVector[39])
                 + z * (d4 * Ay[2] - state.pVector[40])))
            * (-a);
      }

      state.pVector[7] -= (s0 * state.pVector[3]);
      state.pVector[8] -= (s0 * state.pVector[4]);
      state.pVector[9] -= (s0 * state.pVector[5]);
      state.pVector[10] -= (s0 * state.pVector[42]);
      state.pVector[11] -= (s0 * state.pVector[43]);
      state.pVector[12] -= (s0 * state.pVector[44]);

      state.pVector[14] -= (s1 * state.pVector[3]);
      state.pVector[15] -= (s1 * state.pVector[4]);
      state.pVector[16] -= (s1 * state.pVector[5]);
      state.pVector[17] -= (s1 * state.pVector[42]);
      state.pVector[18] -= (s1 * state.pVector[43]);
      state.pVector[19] -= (s1 * state.pVector[44]);

      state.pVector[21] -= (s2 * state.pVector[3]);
      state.pVector[22] -= (s2 * state.pVector[4]);
      state.pVector[23] -= (s2 * state.pVector[5]);
      state.pVector[24] -= (s2 * state.pVector[42]);
      state.pVector[25] -= (s2 * state.pVector[43]);
      state.pVector[26] -= (s2 * state.pVector[44]);

      state.pVector[28] -= (s3 * state.pVector[3]);
      state.pVector[29] -= (s3 * state.pVector[4]);
      state.pVector[30] -= (s3 * state.pVector[5]);
      state.pVector[31] -= (s3 * state.pVector[42]);
      state.pVector[32] -= (s3 * state.pVector[43]);
      state.pVector[33] -= (s3 * state.pVector[44]);

      state.pVector[35] -= (s4 * state.pVector[3]);
      state.pVector[36] -= (s4 * state.pVector[4]);
      state.pVector[37] -= (s4 * state.pVector[5]);
      state.pVector[38] -= (s4 * state.pVector[42]);
      state.pVector[39] -= (s4 * state.pVector[43]);
      state.pVector[40] -= (s4 * state.pVector[44]);

      double P3, P4,
          C = state.pVector[3] * state.pVector[3]
          + state.pVector[4] * state.pVector[4];
      if (C > 1.e-20) {
        C  = 1. / C;
        P3 = state.pVector[3] * C;
        P4 = state.pVector[4] * C;
        C  = -sqrt(C);
      } else {
        C  = -1.e10;
        P3 = 1.;
        P4 = 0.;
      }

      double MA[3] = {Ax[0], Ax[1], Ax[2]};
      double MB[3] = {Ay[0], Ay[1], Ay[2]};
      // Jacobian production of transport and to_local
      if (surface.type() == Surface::Disc) {
        // the vector from the disc surface to the p
        const auto& sfc  = surface.center();
        double      d[3] = {state.pVector[0] - sfc(0),
                       state.pVector[1] - sfc(1),
                       state.pVector[2] - sfc(2)};
        // this needs the transformation to polar coordinates
        double RC = d[0] * Ax[0] + d[1] * Ax[1] + d[2] * Ax[2];
        double RS = d[0] * Ay[0] + d[1] * Ay[1] + d[2] * Ay[2];
        double R2 = RC * RC + RS * RS;

        // inverse radius
        double Ri = 1. / sqrt(R2);
        MA[0]     = (RC * Ax[0] + RS * Ay[0]) * Ri;
        MA[1]     = (RC * Ax[1] + RS * Ay[1]) * Ri;
        MA[2]     = (RC * Ax[2] + RS * Ay[2]) * Ri;
        MB[0] = (RC * Ay[0] - RS * Ax[0]) * (Ri = 1. / R2);
        MB[1] = (RC * Ay[1] - RS * Ax[1]) * Ri;
        MB[2] = (RC * Ay[2] - RS * Ax[2]) * Ri;
      }

      state.jacobian[0] = MA[0] * state.pVector[7] + MA[1] * state.pVector[8]
          + MA[2] * state.pVector[9];  // dL0/dL0
      state.jacobian[1] = MA[0] * state.pVector[14] + MA[1] * state.pVector[15]
          + MA[2] * state.pVector[16];  // dL0/dL1
      state.jacobian[2] = MA[0] * state.pVector[21] + MA[1] * state.pVector[22]
          + MA[2] * state.pVector[23];  // dL0/dPhi
      state.jacobian[3] = MA[0] * state.pVector[28] + MA[1] * state.pVector[29]
          + MA[2] * state.pVector[30];  // dL0/dThe
      state.jacobian[4] = MA[0] * state.pVector[35] + MA[1] * state.pVector[36]
          + MA[2] * state.pVector[37];  // dL0/dCM

      state.jacobian[5] = MB[0] * state.pVector[7] + MB[1] * state.pVector[8]
          + MB[2] * state.pVector[9];  // dL1/dL0
      state.jacobian[6] = MB[0] * state.pVector[14] + MB[1] * state.pVector[15]
          + MB[2] * state.pVector[16];  // dL1/dL1
      state.jacobian[7] = MB[0] * state.pVector[21] + MB[1] * state.pVector[22]
          + MB[2] * state.pVector[23];  // dL1/dPhi
      state.jacobian[8] = MB[0] * state.pVector[28] + MB[1] * state.pVector[29]
          + MB[2] * state.pVector[30];  // dL1/dThe
      state.jacobian[9] = MB[0] * state.pVector[35] + MB[1] * state.pVector[36]
          + MB[2] * state.pVector[37];  // dL1/dCM

      state.jacobian[10]
          = P3 * state.pVector[11] - P4 * state.pVector[10];  // dPhi/dL0
      state.jacobian[11]
          = P3 * state.pVector[18] - P4 * state.pVector[17];  // dPhi/dL1
      state.jacobian[12]
          = P3 * state.pVector[25] - P4 * state.pVector[24];  // dPhi/dPhi
      state.jacobian[13]
          = P3 * state.pVector[32] - P4 * state.pVector[31];  // dPhi/dThe
      state.jacobian[14]
          = P3 * state.pVector[39] - P4 * state.pVector[38];  // dPhi/dCM
      state.jacobian[15] = C * state.pVector[12];             // dThe/dL0
      state.jacobian[16] = C * state.pVector[19];             // dThe/dL1
      state.jacobian[17] = C * state.pVector[26];             // dThe/dPhi
      state.jacobian[18] = C * state.pVector[33];             // dThe/dThe
      state.jacobian[19] = C * state.pVector[40];             // dThe/dCM
      state.jacobian[20] = 0.;                                // dCM /dL0
      state.jacobian[21] = 0.;                                // dCM /dL1
      state.jacobian[22] = 0.;                                // dCM /dPhi
      state.jacobian[23] = 0.;                                // dCM /dTheta
      state.jacobian[24] = state.pVector[41];                 // dCM /dCM

      Eigen::
          Map<Eigen::Matrix<double, NGlobalPars, NGlobalPars, Eigen::RowMajor>>
              J(state.jacobian);

      cov = std::make_unique<const ActsSymMatrixD<NGlobalPars>>(
          J * (*state.covariance) * J.transpose());
    }

    // Fill the end parameters
    result.endParameters = std::make_unique<const BoundParameters>(
        std::move(cov), gp, mom, charge(state), surface.getSharedPtr());
  }

  AtlasSTEPStepper(bfield_t bField = bfield_t()) : m_bField(std::move(bField)){};

  Vector3D
  getField(State& state, const Vector3D& pos) const
  {
    // get the field from the cell
    state.field = m_bField.getField(pos, state.fieldCache);
    return state.field;
  }

  inline double 
  dEdl_ionization(double p, std::shared_ptr<const Material> mat, int particle, double particleMass, double& sigma, double& kazL) const {

    const double factor = (1./3.59524); // the compiler will evaulate this

    double path = 1.; // this is a scaling factor for the landau convolution

    sigma = 0.;
    if ( mat->Z()<1 ) return 0.;    

    double Ionization = 0.;
    
    // kinetic variables
    // and the electron mass in MeV
    double me    = 0.51099891 * units::_MeV;
    double m     = particleMass;
    double mfrac = me/m;
    double E     = sqrt(p*p+m*m);
    double beta  = p/E;
    double gamma = E/m;
    
    //Ionization - Bethe-Bloch
    double I = 16.e-6 * std::pow(mat->Z(),0.9); //16 eV * Z**0.9 - bring to MeV
    
    //K/A*Z = 0.5 * 30.7075MeV/(g/mm2) * Z/A * rho[g/mm3]  / scale to mm by this
    double kaz = 0.5*30.7075*mat->zOverAtimesRho();

    //  sigmaL of Landau 
    sigma = 4*kaz*beta/beta;       // dsigma/dl
    kazL = 0.;
  
    double MOP = 0.;
    if (particle == 11){
      // for electrons use slightly different BetheBloch adaption
      // see Stampfer, et al, "Track Fitting With Energy Loss", Comp. Pyhs. Comm. 79 (1994), 157-164
      Ionization = -kaz*(2.*log(2.*me/I)+3.*log(gamma) - 1.95);
    }  else {

      double eta2 = beta*gamma; eta2 *= eta2;
      // density effect, only valid for high energies (gamma > 10 -> p > 1GeV for muons)
      double delta = 0.;
      if (gamma > 10.) {
	      double eplasma = 28.816e-6 * sqrt(1000.*mat->zOverAtimesRho());
	      delta = 2.*log(eplasma/I) + log(eta2) - 1.;
      }
      // tmax - cut off energy
      double tMax = 2.*eta2*me/(1.+2.*gamma*mfrac+mfrac*mfrac);
      // divide by beta^2 for non-electrons
      kaz /= beta*beta;
      Ionization = -kaz*(log(2.*me*eta2*tMax/(I*I)) - 2.*(beta*beta) - delta);
      //
      // the landau sigmaL is path length dependent - take what's given
      //
      //    PDG formula 27.11 for MOP value from http://pdg.lbl.gov/2011/reviews/rpp2011-rev-passage-particles-matter.pdf 
      //
      MOP =  -kaz*(log(2.*me*eta2/I) + log(path*kaz/I) + 0.2 - (beta*beta) - delta);
      sigma = -(Ionization - MOP)*factor; 
      kazL = kaz*factor;
    }
    // return mean or mop
    return Ionization;
  }

  inline double 
  dEdl_radiation(double p, std::shared_ptr<const Material> mat, int particle, double particleMass, double& sigma) const{
    sigma = 0.;
    if ( mat->X0()==0. ) return 0.;    

    // preparation of kinetic constants
    double m     = particleMass;
    double me    = 0.51099891 * units::_MeV;
    double mfrac = me/m;
    double E     = sqrt(p*p+m*m);

    //Bremsstrahlung - Bethe-Heitler
    double Radiation = -E*mfrac*mfrac;
    // sigma is rms of steep exponential part of radiation 
    sigma = -Radiation;

    if ((particle == 13) && (E > 8000.)) {
      if (E < 1.e6) {
    	      Radiation += 0.5345 - 6.803e-5*E - 2.278e-11*E*E + 9.899e-18*E*E*E; //E below 1 TeV
          sigma     += (0.1828 -3.966e-3*sqrt(E) + 2.151e-5*E ); // idem
      } else {
    	      Radiation += 2.986 - 9.253e-5*E; //E above 1 TeV
          sigma     += 17.73 + 2.409e-5*(E-1000000.); // idem
      }
    }

    sigma = sigma/mat->X0();  
    
    return Radiation/mat->X0();  
  }
  
template <typename propagator_state_t>
double 
dEds(propagator_state_t& state, double p) const
{
  state.stepping.m_delIoni = 0.; 
  state.stepping.m_delRad = 0.; 
  state.stepping.m_kazL = 0.;

  if (state.navigation.currentVolume == nullptr || state.navigation.currentVolume->material()->X0() == 0 || state.navigation.currentVolume->material()->Z() ==0) return 0.; 

  state.stepping.m_delIoni = dEdl_ionization(p, state.navigation.currentVolume->material(), state.options.absPdgCode, state.options.mass, state.stepping.m_sigmaIoni, state.stepping.m_kazL);    

  state.stepping.m_delRad  = dEdl_radiation(p, state.navigation.currentVolume->material(), state.options.absPdgCode, state.options.mass, state.stepping.m_sigmaRad);

  double eLoss = state.stepping.m_MPV ? 0.9*state.stepping.m_delIoni + 0.15*state.stepping.m_delRad : state.stepping.m_delIoni + state.stepping.m_delRad;

  return eLoss;
}  
  
  /// Perform the actual step on the state
  ///
  /// @param state is the provided stepper state (caller keeps thread locality)
  template <typename propagator_state_t>
  double
  step(propagator_state_t& state) const
  {
  double sol = 0.0299792458;			// Speed of light
  double charge;
  state.stepping.pVector[6] >= 0. ? charge = 1. : charge = -1.;      // Set charge
  double     lambda1 = state.stepping.pVector[6], lambda2 = state.stepping.pVector[6];	// Store inverse momentum for Jacobian transport
  double     lambda3 = state.stepping.pVector[6], lambda4 = state.stepping.pVector[6];
  double     dP1=0., dP2=0., dP3=0., dP4=0.;    // dp/ds = -g*E/p for positive g=dE/ds
  double     dL1=0., dL2=0., dL3=0., dL4=0.;    // factor used for calculating dCM/dCM, pVector[41], in the Jacobian.
  double     initialMomentum = fabs( 1./state.stepping.pVector[6]);  // Set initial momentum
  Vector3D initialPos( state.stepping.pVector[0], state.stepping.pVector[1], state.stepping.pVector[2]);	// Set initial values for position
  Vector3D initialDir( state.stepping.pVector[3], state.stepping.pVector[4], state.stepping.pVector[5]);	// Set initial values for direction.
  Vector3D dir1, dir2, dir3, dir4;            // Directions at the different points. Used by the error propagation
  ActsVectorF<12> BG23 = ActsVectorF<12>::Zero();
  ActsVectorF<12> BG4 = ActsVectorF<12>::Zero();
  double     g = 0.;                            // Energyloss in Mev/mm.
  double     dgdl = 0.;                         // dg/dlambda in Mev/mm.
  double     particleMass = state.options.mass; //Get particle mass from ParticleHypothesis
  int        steps = 0;
  double h = state.stepping.stepSize;

  //POINT 1. Only calculate this once per step, even if step is rejected. This point is independant of the step length, h
  double     momentum = initialMomentum;        // Current momentum
  if (state.stepping.m_energyLoss && state.navigation.currentVolume != nullptr && state.navigation.currentVolume->material() != nullptr) {
    g = dEds( state, momentum); //Use the same energy loss throughout the step.
    double E = std::sqrt( momentum*momentum + particleMass*particleMass);
    dP1 = g*E/momentum;
    if (state.stepping.errorPropagation) {
		      if (state.stepping.m_includeGgradient) {} // TODO
      dL1 = -lambda1*lambda1*g*E*(3.-(momentum*momentum)/(E*E)) - lambda1*lambda1*lambda1*E*dgdl;
    }
  }

  Vector3D pos = initialPos; // Coordinate solution
  Vector3D dir = initialDir; // Direction solution
  dir1 = dir;
  if (state.stepping.firstStep) {	       // Poll BG1 if this is the first step, else use recycled BG4
    state.stepping.firstStep = false;
    state.stepping.field = m_bField.getField(pos, state.stepping.fieldCache);
    state.stepping.BG1 << state.stepping.field.x(), state.stepping.field.y(), state.stepping.field.z(),0.,0.,0.,0.,0.,0.,0.,0.,0.; //Get the gradients needed for the error propagation if errorPropagation=true
  }
  // Lorentz force, d2r/ds2 = lambda * (dir x B)
  Vector3D k1( dir.y()*state.stepping.BG1[2] - dir.z()*state.stepping.BG1[1], dir.z()*state.stepping.BG1[0] - dir.x()*state.stepping.BG1[2], dir.x()*state.stepping.BG1[1] - dir.y()*state.stepping.BG1[0]);
  k1 = sol * lambda1 * k1;
  while (true) { //Repeat step until error estimate is within the requested tolerance
std::cout << steps << std::endl;
    if ((unsigned int) steps++ > state.options.maxSteps) return 0.; //Abort propagation
    //POINT 2
    if (state.stepping.m_energyLoss && state.navigation.currentVolume != nullptr && state.navigation.currentVolume->material() != nullptr) {
      momentum = initialMomentum + (h/2.)*dP1;
      if (momentum <= state.stepping.m_momentumCutOff) { h *= 0.5; state.stepping.stepSize.update(h, detail::ConstrainedStep::accuracy); continue;} //Abort propagation
      double E = std::sqrt( momentum*momentum + particleMass*particleMass);
      dP2 = g*E/momentum;
      lambda2 = charge/momentum;
      if (state.stepping.errorPropagation) {
        dL2 = -lambda2*lambda2*g*E*(3.-(momentum*momentum)/(E*E)) - lambda2*lambda2*lambda2*E*dgdl;
      }
    }

    pos = initialPos + (h/2.)*initialDir + (h*h/8.)*k1;
    dir = initialDir + (h/2.)*k1;
    dir2 = dir;
    state.stepping.field = m_bField.getField(pos, state.stepping.fieldCache);
    BG23 << state.stepping.field.x(), state.stepping.field.y(), state.stepping.field.z(),0.,0.,0.,0.,0.,0.,0.,0.,0.;
    Vector3D k2( dir.y()*BG23[2] - dir.z()*BG23[1], dir.z()*BG23[0] - dir.x()*BG23[2],
      dir.x()*BG23[1] - dir.y()*BG23[0]);
    k2 = sol * lambda2 * k2;
   
    //POINT 3. Same position and B-field as point 2.
    if (state.stepping.m_energyLoss && state.navigation.currentVolume != nullptr && state.navigation.currentVolume->material() != nullptr) {
      momentum = initialMomentum + (h/2.)*dP2;
      if (momentum <= state.stepping.m_momentumCutOff) { h *= 0.5; state.stepping.stepSize.update(h, detail::ConstrainedStep::accuracy); continue;} //Abort propagation
      double E = std::sqrt( momentum*momentum + particleMass*particleMass);
      dP3 = g*E/momentum;
      lambda3 = charge/momentum;
      if (state.stepping.errorPropagation) {
	dL3 = -lambda3*lambda3*g*E*(3.-(momentum*momentum)/(E*E)) - lambda3*lambda3*lambda3*E*dgdl;
      }
    }
    dir = initialDir + (h/2.)*k2;
    dir3 = dir;
    Vector3D k3( dir.y()*BG23[2] - dir.z()*BG23[1], dir.z()*BG23[0] - dir.x()*BG23[2],
      dir.x()*BG23[1] - dir.y()*BG23[0]);
    k3 = sol * lambda3 * k3;

    //POINT 4
    if (state.stepping.m_energyLoss && state.navigation.currentVolume != nullptr && state.navigation.currentVolume->material() != nullptr) {
      momentum = initialMomentum + h*dP3;
      if (momentum <= state.stepping.m_momentumCutOff) { h *= 0.5; state.stepping.stepSize.update(h, detail::ConstrainedStep::accuracy); continue;} //Abort propagation
      double E = std::sqrt( momentum*momentum + particleMass*particleMass);
      dP4 = g*E/momentum;
      lambda4 = charge/momentum;
      if (state.stepping.errorPropagation) {
	dL4 = -lambda4*lambda4*g*E*(3.-(momentum*momentum)/(E*E)) - lambda4*lambda4*lambda4*E*dgdl;
      }
    }
    pos = initialPos + h*initialDir + (h*h/2.)*k3;
    dir = initialDir + h*k3;
    dir4 = dir;
    state.stepping.field = m_bField.getField(pos, state.stepping.fieldCache);
    BG4 << state.stepping.field.x(), state.stepping.field.y(), state.stepping.field.z(),0.,0.,0.,0.,0.,0.,0.,0.,0.;
    Vector3D k4( dir.y()*BG4[2] - dir.z()*BG4[1], dir.z()*BG4[0] - dir.x()*BG4[2], dir.x()*BG4[1] - dir.y()*BG4[0]);
    k4 = sol * lambda4 * k4;

    //Estimate local error and avoid division by zero
    Vector3D errorPos((h*h) * (k1 - k2 - k3 + k4));
    double errorEstimate = std::max( errorPos.norm(), 1e-20);

    //Use the error estimate to calculate new step length. h is returned by reference.
    double distanceStepped = h; //Store old step length.
    h = h*std::min( std::max( 0.25, std::pow((state.options.tolerance / errorEstimate), 0.25)), 4.);

    //Repeat step with the new step size if error is too big.
    if (errorEstimate > 4.*state.options.tolerance) continue;

    //Update inverse momentum if energyloss is switched on
    if (state.stepping.m_energyLoss && state.navigation.currentVolume != nullptr && state.navigation.currentVolume->material() != nullptr) {
      if (initialMomentum + (distanceStepped/6.)*(dP1 + 2.*dP2 + 2.*dP3 + dP4) <= state.stepping.m_momentumCutOff) { h *= 0.5; state.stepping.stepSize.update(h, detail::ConstrainedStep::accuracy); continue;} //Abort propagation
      momentum = initialMomentum + (distanceStepped/6.)*(dP1 + 2.*dP2 + 2.*dP3 + dP4);
      state.stepping.pVector[6] = charge/momentum;
    }
    
    //Step was ok. Store solutions.
    //Update positions.
    pos = initialPos + distanceStepped*initialDir + (distanceStepped*distanceStepped/6.)*(k1 + k2 + k3);

    state.stepping.pVector[0] = pos.x();
    state.stepping.pVector[1] = pos.y();
    state.stepping.pVector[2] = pos.z();

    //update directions
    dir = initialDir + (distanceStepped/6.)*(k1 + 2.*k2 + 2.*k3 + k4);
    state.stepping.pVector[3] = dir.x();
    state.stepping.pVector[4] = dir.y();
    state.stepping.pVector[5] = dir.z();

    //Normalize direction
    double norm = 1./std::sqrt( state.stepping.pVector[3]*state.stepping.pVector[3] + state.stepping.pVector[4]*state.stepping.pVector[4] + state.stepping.pVector[5]*state.stepping.pVector[5]);
    state.stepping.pVector[3] = norm*state.stepping.pVector[3];
    state.stepping.pVector[4] = norm*state.stepping.pVector[4];
    state.stepping.pVector[5] = norm*state.stepping.pVector[5];




    //dDir provides a small correction to the final tiny step in PropagateWithJacobian
    Vector3D dDir;
    dDir[0] = k4.x();
    dDir[1] = k4.y();
    dDir[2] = k4.z();

    //Transport Jacobian using the same step length, points and magnetic fields as for the track parameters
    //BG contains Bx, By, Bz, dBx/dx, dBx/dy, dBx/dz, dBy/dx, dBy/dy, dBy/dz, dBz/dx, dBz/dy, dBz/dz
    //               0   1   2   3       4       5       6       7       8       9       10      11
    if (state.stepping.errorPropagation) {
      double     initialjLambda = state.stepping.pVector[41]; //dLambda/dLambda
      double     jLambda = initialjLambda;
      double     jdL1=0., jdL2=0., jdL3=0., jdL4=0.;

      for (int i=7; i<42; i+=7) {

        //POINT 1
	Vector3D initialjPos( state.stepping.pVector[i], state.stepping.pVector[i+1], state.stepping.pVector[i+2]);
	Vector3D jPos = initialjPos;
  	Vector3D initialjDir( state.stepping.pVector[i+3], state.stepping.pVector[i+4], state.stepping.pVector[i+5]);
  	Vector3D jDir = initialjDir;

	//B-field terms
    Vector3D jk1( jDir.y()*state.stepping.BG1[2] - jDir.z()*state.stepping.BG1[1], jDir.z()*state.stepping.BG1[0] - jDir.x()*state.stepping.BG1[2],
	  jDir.x()*state.stepping.BG1[1] - jDir.y()*state.stepping.BG1[0]);

        //B-field gradient terms
	if (state.stepping.m_includeBgradients) {}
        jk1 = sol * lambda1 * jk1;

	//Last column of the A-matrix
        if (i==35) {                  //Only dLambda/dLambda, in the last row of the jacobian, is different from zero
          jdL1 = dL1*jLambda;         //Energy loss term. dL1 = 0 if no material effects -> jLambda = pVector[41] is unchanged
	  jk1 = jk1 + k1*jLambda/lambda1; //B-field terms
	}

	//POINT 2
	jPos = initialjPos + (distanceStepped/2.)*initialjDir + (distanceStepped*distanceStepped/8.)*jk1;
	jDir = initialjDir + (distanceStepped/2.)*jk1;

       Vector3D jk2( jDir.y()*BG23[2] - jDir.z()*BG23[1], jDir.z()*BG23[0] - jDir.x()*BG23[2],
	  jDir.x()*BG23[1] - jDir.y()*BG23[0]);

	if (state.stepping.m_includeBgradients) {}
    	jk2 = sol * lambda2 * jk2;

	if (i==35) {
	  jLambda = initialjLambda + (distanceStepped/2.)*jdL1;
          jdL2 = dL2*jLambda;
	  jk2 = jk2 + k2*jLambda/lambda2;
	}

	//POINT 3
    	jDir = initialjDir + (distanceStepped/2.)*jk2;

    	Vector3D jk3( jDir.y()*BG23[2] - jDir.z()*BG23[1], jDir.z()*BG23[0] - jDir.x()*BG23[2],
	  jDir.x()*BG23[1] - jDir.y()*BG23[0]);

	if (state.stepping.m_includeBgradients) {}
    	jk3 = sol * lambda3 * jk3;

	if (i==35) {
	  jLambda = initialjLambda + (distanceStepped/2.)*jdL2;
          jdL3 = dL3*jLambda;
	  jk3 = jk3 + k3*jLambda/lambda3;
        }

    	//POINT 4
	jPos = initialjPos + distanceStepped*initialjDir + (distanceStepped*distanceStepped/2.)*jk3;
    	jDir = initialjDir + distanceStepped*jk3;

     	Vector3D jk4( jDir.y()*BG4[2] - jDir.z()*BG4[1], jDir.z()*BG4[0] - jDir.x()*BG4[2],
	  jDir.x()*BG4[1] - jDir.y()*BG4[0]);

	if (state.stepping.m_includeBgradients) {}
    	jk4 = sol * lambda4 * jk4;

        if (i==35) {
	  jLambda = initialjLambda + distanceStepped*jdL3;
          jdL4 = dL4*jLambda;
	  jk4 = jk4 + k4*jLambda/lambda4;
        }

        //solution
	jPos = initialjPos + distanceStepped*initialjDir + (distanceStepped*distanceStepped/6.)*(jk1 + jk2 + jk3);
        jDir = initialjDir + (distanceStepped/6.)*(jk1 + 2.*jk2 + 2.*jk3 + jk4);
	if (i==35) {
	  jLambda = initialjLambda + (distanceStepped/6.)*(jdL1 + 2.*jdL2 + 2.*jdL3 + jdL4);
        }

        //Update positions
        state.stepping.pVector[i]   = jPos.x();
        state.stepping.pVector[i+1] = jPos.y();
        state.stepping.pVector[i+2] = jPos.z();

        //update directions
        state.stepping.pVector[i+3] = jDir.x();
        state.stepping.pVector[i+4] = jDir.y();
        state.stepping.pVector[i+5] = jDir.z();
      }
      state.stepping.pVector[41] = jLambda; //update dLambda/dLambda
    }

    //Store BG4 for use as BG1 in the next step
    for (int i=0; i<12; i++) {
      state.stepping.BG1[i] = BG4[i];
    }
    return h;
  }
  }

private:
  bfield_t m_bField;
};

}  // namespace Acts