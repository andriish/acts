// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"
#include <random>
#include "ActsFatras/EventData/ProcessType.hpp"
#include <vector>
#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {   
  class PhotonConversion {
  public:
  
      /** The cut from which on the child products are followed */
      double                                       m_minChildEnergy = 50. * Acts::UnitConstants::MeV;
      double                                       m_childEnergyScaleFactor = 2.;
      double                                       m_conversionProbScaleFactor = 0.98;
      
      /** interface for processing of the pair production */
      bool pairProduction(const Trk::MaterialProperties& mprop,
			  double pathCorrection,
			  double p) const;
      
      /** interface for processing of the presampled pair production */      
      /** interface for processing of the presampled nuclear interactions on layer*/
      template <typename generator_t>
      std::vector<Particle> doConversion(generator_t& generator ,double time, const Particle& parm) const;
    
   private:
      /** record childs - create interface already for mu-/mu+ pair production*/
      template <typename generator_t>
      std::vector<Particle> recordChilds(generator_t& generator, const Particle& photon, 
			double childEnergy,
			const Acts::Vector3D& childDirection) const;
	
      /** simulate the child energy */
      template <typename generator_t>
      Particle::Scalar childEnergyFraction(generator_t& generator, Particle::Scalar gammaMom) const;
      
      /** simulate the one child direction - the second one is given clear then*/
      template <typename generator_t>
      Particle::Vector3 childDirection(generator_t& generator, const Particle::Vector4& gammaMom4) const;
      
      /** helper functions for the Phi1/phi2 */
      double phi1(double delta) const;
      
      /** helper functions for the Phi1/phi2 */
      double phi2(double delta) const;
            
      /** Inverse fine structure constant */
      constexpr double                                s_alphaEM = 1. / 137.;
      constexpr double s_alphaEM2 = s_alphaEM * s_alphaEM;
      constexpr double                                s_oneOverThree = 1. / 3.;
   };
  

inline double PhotonConversionTool::phi1(double delta) const {
  if (delta <= 1.)
     return 20.867 - 3.242 * delta  + 0.625 * delta * delta;
  return 21.12 - 4.184*log(delta+0.952);
}

inline double PhotonConversionTool::phi2(double delta) const {
  if (delta <= 1.)
     return 20.209 - 1.930 * delta  - 0.086*delta*delta;
   return 21.12 - 4.184*log(delta+0.952);
}

}


#include <math.h>

template <typename generator_t>
Particle::Scalar ActsFatras::PhotonConversionTool::pairProduction(generator_t& generator, Particle::Scalar momentum) const
{

 // use for the moment only Al data - Yung Tsai - Rev.Mod.Particle Physics Vol. 46, No.4, October 1974
 // optainef from a fit given in the momentum range 100 10 6 2 1 0.6 0.4 0.2 0.1 GeV

 //// Quadratic background function
 //  Double_t fitFunction(Double_t *x, Double_t *par) {
 //  return par[0] + par[1]*pow(x[0],par[2]);
 // }
 
 // EXT PARAMETER                                   STEP         FIRST
 // NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
 //  1  p0          -7.01612e-03   8.43478e-01   1.62766e-04   1.11914e-05
 //  2  p1           7.69040e-02   1.00059e+00   8.90718e-05  -8.41167e-07
 //  3  p2          -6.07682e-01   5.13256e+00   6.07228e-04  -9.44448e-07

 constexpr double  p0  =       -7.01612e-03;
 constexpr double  p1  =        7.69040e-02;
 constexpr double  p2  =       -6.07682e-01;
 // calculate xi
 const double xi = p0 + p1*pow(p,p2);
 // now calculate what's left
 //~ double attenuation = exp( -7.777e-01*pathCorrection*mprop.thicknessInX0()*(1.-xi) ); // eq. 3.75

  //~ return (m_conversionProbScaleFactor*CLHEP::RandFlat::shoot(m_randomEngine) > attenuation) ? true : false;
  
  //~ return (CLHEP::RandFlat::shoot(m_randomEngine) < 1 - attenuation / m_conversionProbScaleFactor) ? true : false;
  //~ m_conversionProbScaleFactor(rnd - 1) < -attenuation
  //~ m_conversionProbScaleFactor(1 - rnd) > attenuation
  //~ ln(m_conversionProbScaleFactor(1 - rnd)) > -7/9*pathCorrection*mprop.thicknessInX0()*(1.-xi)
  
  constexpr double NineOverSeven = 9 / 7;
  std::uniform_real_distribution<double> uniformDistribution {0., 1.};
  return -NineOverSeven * ln(m_conversionProbScaleFactor * (1 - uniformDistribution(generator))) / (1.-xi);
}

template <typename generator_t>
ActsFatras::Particle::Scalar ActsFatras::PhotonConversionTool::childEnergyFraction(generator_t& generator, Particle::Scalar gammaMom) const {

  /// Some (if not all) calculations are from G4PairProductionRelModel
  
  // the fraction
  const double epsilon0      = findMass(eElectron) / gammaMom;
  // some needed manipolations
  // TODO: why fixed Z? and why is it Al?
  constexpr double Z             = 13.; //mprop.averageZ(); // TODO: the other part was never used
  constexpr double oneOverZpow   = 1. / pow(Z,s_oneOverThree);
  constexpr double alphaZsquare  = s_alphaEM2 * Z * Z;
  // now f(Z) - from Z and s_alpha
  constexpr double fZ            = alphaZsquare * (1. / (1. + alphaZsquare) + 0.20206 - 0.0369 * alphaZsquare + 0.0083 * alphaZsquare * alphaZsquare);
  // delta_max
  constexpr double deltaMax      = exp((42.24 - fZ) * 0.1195) - 0.952;
  // delta_min
  const double deltaMin          = 4. * epsilon0 * 136. * oneOverZpow; 
  // the minimum fraction
  const double epsilon1          = 0.5 - 0.5 * sqrt(1. - deltaMin / deltaMax);
  const double epsilonMin    = std::max(epsilon1, epsilon0);
  // calculate phi1 / phi2 - calculate from deltaMin
  const double Phi1          = phi1(deltaMin);
  const double Phi2          = phi2(deltaMin);
  // then calculate F10/F20
  const double F10           = 3. * Phi1 - Phi2 - fZ;
  const double F20           = 1.5 * Phi1 - 0.5 * Phi2 - fZ;
  // and finally calucate N1, N2
  const double N1            = (0.25 - epsilonMin + epsilonMin * epsilonMin) * F10;
  const double N2            = 1.5 * F20;
  
  // ------------ decide wich one to take
  double epsilon;
  double delta;
  double F;
  const double deltaNumerator = onOverZpow * 136. * epsilon0;
  std::uniform_real_distribution<double> uniformDistribution {0., 1.};
  if ( N1 < uniformDistribution(generator) * (N1 + N2) ) {
    // sample from f1,g1 distribution
    do{
      epsilon = 0.5 - (0.5 - epsilonMin) * pow(uniformDistribution(generator),s_oneOverThree);
      // prepare for the rejection check
      delta   = deltaNumerator / (epsilon - epsilon * epsilon);
      F = 3. * phi1(delta) - phi2(delta) - fZ;   
      // reject ? - or redo the exercise 
    } while(F <= uniformDistribution(generator) * F10);
  } else {
    // sample from f2,g2 distribution
    do {
      epsilon = epsilonMin + (0.5-epsilonMin)*uniformDistribution(generator);
      // prepare for the rejection check
      delta   = deltaNumerator / (epsilon - epsilon * epsilon);
      F = 1.5 * phi1(delta) - 0.5 * phi2(delta) - fZ;   
     // reject ? - or redo the exercise 
    } while (F <= uniformDistribution(generator) * F20);    
  }
  return m_childEnergyScaleFactor * epsilon;
}

template <typename generator_t>
Particle::Vector3 ActsFatras::PhotonConversionTool::childDirection(generator_t& generator, const ActsFatras::Particle::Vector4& gammaMom4) const
{
    // --------------------------------------------------
    // Following the Geant4 approximation from L. Urban
    // the azimutal angle
    
    // the start of the equation
    Particle::Scalar theta = findMass(eElectron) / gammaMom4[eEnergy];

	std::unitform_real_distribution<Particle::Scalar> uniformDistribution {0., 1.};
    const Particle::Scalar u =  -log(uniformDistribution(generator) * uniformDistribution(generator)) * 1.6;
    
    theta *= (uniformDistribution(generator) < 0.25 ) ? u : u*s_oneOverThree; // 9./(9.+27) = 0.25

    // more complex but "more true"
    const Particle::Vector3 gammaMomHep = gammaMom4.template segment<3>(eDir0);
    const Particle::Vector3 newDirectionHep(gammaMomHep.normalized());
    
    // if it runs along the z axis - no good ==> take the x axis
    const Particle::Scalar x = (newDirectionHep[eZ] * newDirectionHep[eZ] > 0.999999) ? 1. : -newDirectionHep[eY];
    const Particle::Scalar y = newDirectionHep[eX];
    
    // deflector direction
    const Particle::Vector3 deflectorHep(x, y, 0.);
    // rotate the new direction for scattering
    Eigen::Transform rotTheta = Eigen::AngleAxis<Particle::Scalar>(theta, deflectorHep);
    
    // and arbitrarily in psi
    Eigen::Transform rotPsi = Eigen::AngleAxis<Particle::Scalar>(uniformDistribution(generator) * 2. * M_PI, gammaMomHep);

    return rotPsi * rotTheta * newDirectionHep;
}

template <typename generator_t>
void ActsFatras::PhotonConversionTool::recordChilds(generator_t& generator, ActsFatras::Particle& photon,
                                                ActsFatras::Particle::Scalar childEnergy,
                                                const ActsFatras::Particle::Vector3& childDirection,
                                                Acts::PdgParticle pdgProduced) const
{
    // Calculate the child momentum
    const float massChild = findMass(pdgProduced);
    const Particle::Scalar momentum1 = sqrt(childEnergy * childEnergy - massChild * massChild);    

    // now properly : energy-momentum conservation
    const Particle::Vector3 vtmp = photom.momentum4().template segment<3>(Acts::eDir0) - momentum1 * childDirection;
    const Particle::Scalar momentum2 = vtmp.norm();

    // charge sampling
    std::uniform_int_distribution<> uniformDistribution {0, 1};
    Particle::Scalar charge1;
    Particle::Scalar charge2;
    if (uniformDistribution(generator) == 1) {
      charge1 = -1.;
      charge2 =  1.;
    } else {
      charge1 =  1.;
      charge2 = -1.;
    }

    // add the new secondary states to the ISF particle stack
    const Acts::PdgParticle    pdg1  = std::copysign(Acts::PdgParticle::eElectron, charge1);
    const Acts::PdgParticle    pdg2  = std::copysign(Acts::PdgParticle::eElectron, charge2);

    // remove soft children
    int nchild = 0;
    if ( p1 > m_minChildEnergy ) nchild++;
    if ( p2 > m_minChildEnergy ) nchild++;

    std::vector<Particle> children;
    children.reserve(2);

    if (  momentum1 > m_minChildEnergy ) {
		Particle child = Particle(Barcode(), pdg1).setPosition4(photon.position4()).setDirection(childDirection).setAbsMomentum(p1).setProcess(ProcessType::ePhotonConversion);
	  children.push_back(std::move(child));
    }

    if (  momentum2 > m_minChildEnergy ) {
		Particle child = Particle(Barcode(), pdg2).setPosition4(photon.position4()).setDirection(childDirection)
			.setAbsMomentum(p2).setProcess(ProcessType::ePhotonConversion);
	  children.push_back(std::move(child));
    }
    return children;
}

template <typename generator_t>
std::vector<ActsFatras::Particle> ActsFatras::PhotonConversionTool::doConversion(generator_t& generator,
	double time, const ActsFatras::Particle& parm) const {
  const double p = parm.absMomentum();

  // get the energy
  const Particle::Scalar childEnergy = p * childEnergyFraction(generator, p);

  // now get the deflection
  Particle::Vector3 childDir = childDirection(generator, parm.momentum4());
  // verbose output
   return recordChilds(parm,
       childEnergy,
       childDir);
}