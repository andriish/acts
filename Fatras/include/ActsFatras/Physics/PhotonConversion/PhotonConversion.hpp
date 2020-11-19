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
#include <math.h>

namespace ActsFatras {
	/// @brief This class handles the photon conversion. It evaluates the distance after which the interaction will occur and the final state due the interaction itself.
  class PhotonConversion {
  public:
  
      /// Lower energy threshold for produced children
      double                                       m_minChildEnergy = 50. * Acts::UnitConstants::MeV;
      /// Scaling factor of children energy
      double                                       m_childEnergyScaleFactor = 2.;
      /// Scaling factor for photon conversion probability
      double                                       m_conversionProbScaleFactor = 0.98;
      
      /// @brief Method for evaluating the distance after which the photon conversion will occur
      ///
      /// @tparam generator_t Type of the random number generator
      /// @param [in, out] generator The random number generator
      /// @param [in] momentum The momentum of the particle
      ///
      /// @param The distance in X_0
		template <typename generator_t>
		Particle::Scalar pairProduction(generator_t& generator, Particle::Scalar momentum) const;
      
      /// @brief This method evaluates the final state due to the photon conversion
      ///
      /// @tparam generator_t Type of the random number generator
      /// @param [in, out] generator The random number generator
      /// @param [in] particle The interacting photon
      ///
      /// @return Vector containing the final state leptons
      template <typename generator_t>
      std::vector<Particle> doConversion(generator_t& generator, const Particle& particle) const;
    
   private:
      /// @brief This method constructs and returns the child particles
      ///
      /// @tparam generator_t Type of the random number generator
      /// @param [in, out] generator The random number generator
      /// @param [in] photon The interacting photon
      /// @param [in] childEnergy The energy of one child particle
      /// @param [in] childDirection The direction of the child particle
      ///
      /// @return Vector containing the produced leptons
      template <typename generator_t>
      std::vector<Particle> recordChilds(generator_t& generator, const Particle& photon, 
			double childEnergy,
			const Acts::Vector3D& childDirection) const;
	
      /// @brief This method evaluates the energy of a child particle
      ///
      /// @tparam generator_t Type of the random number generator
      /// @param [in, out] generator The random number generator
      /// @param [in] gammaMom The momentum of the photon
      ///
      /// @return The energy of the child particle
      template <typename generator_t>
      Particle::Scalar childEnergyFraction(generator_t& generator, Particle::Scalar gammaMom) const;
      
      template <typename generator_t>
      Particle::Scalar childEnergyFractionReDo(generator_t& generator, Particle::Scalar gammaMom) const;
      
      /// @brief This method evaluates the direction of a child particle
      ///
      /// @tparam generator_t Type of the random number generator
      /// @param [in, out] generator The random number generator
      /// @param [in] gammaMom4 The momentum four vector of the photon
      ///
      /// @return The direction vector of the child particle
      template <typename generator_t>
      Particle::Vector3 childDirection(generator_t& generator, const Particle::Vector4& gammaMom4) const;
      
      /// Helper methods for angular evaluation
      double phi1(double delta) const;
      double phi2(double delta) const;
            
      /// Fine structure constant
      constexpr double                                alphaEM = 1. / 137.;
      constexpr double alphaEM2 = alphaEM * alphaEM;
      
      /// Irrational numbers
      constexpr double                                oneOverThree = 1. / 3.;
      constexpr double nineOverSeven = 9 / 7;
      
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
   
    //From G4Element
    constexpr double computeCoulombFactor()
	{
	  //
	  //  Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

	  constexpr double k1 = 0.0083;
	  constexpr double k2 = 0.20206;
	  constexpr double k3 = 0.0020; // This term is missing in Athena
	  constexpr double k4 = 0.0369;

		constexpr float z = 13.; // Aluminium

	  constexpr double az2 = (alphaEM * z)*(alphaEM * z);
	  constexpr double az4 = az2 * az2;

	  return (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
	}

// Compute the value of the screening function 3*PHI1(delta) - PHI2(delta):
inline double ScreenFunction1(const double delta)
{
  return (delta > 1.4) ? 42.038 - 8.29 * log(delta + 0.958) 
                       : 42.184 - delta * (7.444 - 1.623 * delta);
}


// Compute the value of the screening function 1.5*PHI1(delta) +0.5*PHI2(delta):
inline double ScreenFunction2(const double delta)
{
  return (delta > 1.4) ? 42.038 - 8.29 * log(delta + 0.958)
                       : 41.326 - delta * (5.848 - 0.902 * delta);
}
   
   };
  

inline double PhotonConversionTool::phi1(double delta) const {
  /// G4PairProductionRelModel::Phi1
  if (delta <= 1.)
     return 20.867 - 3.242 * delta  + 0.625 * delta * delta;
  return 21.12 - 4.184*log(delta+0.952);
}

inline double PhotonConversionTool::phi2(double delta) const {
  /// G4PairProductionRelModel::Phi2
  if (delta <= 1.)
     return 20.209 - 1.930 * delta  - 0.086*delta*delta;
   return 21.12 - 4.184*log(delta+0.952);
}

}

template <typename generator_t>
Particle::Scalar 
ActsFatras::PhotonConversion::pairProduction(generator_t& generator, Particle::Scalar momentum) const
{
	// eq. 3.75
	
	
 // calculate xi
 const double xi = p0 + p1*pow(p,p2);
 // now calculate what's left
 //~ double attenuation = exp( -7.777e-01*pathCorrection*mprop.thicknessInX0()*(1.-xi) ); 

  //~ return (m_conversionProbScaleFactor*CLHEP::RandFlat::shoot(m_randomEngine) > attenuation) ? true : false;
  
  //~ return (CLHEP::RandFlat::shoot(m_randomEngine) < 1 - attenuation / m_conversionProbScaleFactor) ? true : false;
  //~ m_conversionProbScaleFactor(rnd - 1) < -attenuation
  //~ m_conversionProbScaleFactor(1 - rnd) > attenuation
  //~ ln(m_conversionProbScaleFactor(1 - rnd)) > -7/9*pathCorrection*mprop.thicknessInX0()*(1.-xi)
  
  std::uniform_real_distribution<double> uniformDistribution {0., 1.};
  return -nineOverSeven * ln(m_conversionProbScaleFactor * (1 - uniformDistribution(generator))) / (1.-xi);
}

template <typename generator_t>
ActsFatras::Particle::Scalar 
ActsFatras::PhotonConversion::childEnergyFraction(generator_t& generator, Particle::Scalar gammaMom) const {

  /// Some (if not all) calculations are from G4PairProductionRelModel
  
  // the fraction
  const double epsilon0      = findMass(eElectron) / gammaMom; // == eps0
  // some needed manipolations
  // TODO: why fixed Z? and why is it Al?
  constexpr double Z             = 13.; //mprop.averageZ(); // TODO: the other part was never used
  constexpr double oneOverZpow   = 1. / pow(Z,s_oneOverThree); // == part of gElementData[iZet]->fDeltaFactor
  constexpr double alphaZsquare  = alphaEM2 * Z * Z;
  // now f(Z) - from Z and s_alpha
  constexpr double fZ            = alphaZsquare * (1. / (1. + alphaZsquare) + 0.20206 - 0.0369 * alphaZsquare + 0.0083 * alphaZsquare * alphaZsquare);
  // delta_max
  constexpr double deltaMax      = exp((42.24 - fZ) * 0.1195) + 0.952; //From G4: "F.Hariri : correct sign of last term"
  // delta_min 			 G4Exp((42.038 - FZHigh)/8.29) - 0.958; 	FZHigh = 8.*(logZ13 + fc);
  const double deltaMin          = 4. * epsilon0 * 136. * oneOverZpow; // == deltaMin
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
Particle::Scalar childEnergyFractionReDo(generator_t& generator, Particle::Scalar gammaMom) const {
	
  // Require photon energy >= 100 MeV
  if(gammaMom < 100. * Acts::UnitConstants::MeV)
	return 0.; // TODO: The return value should be different
	
  const double    eps0        = findMass(eElectron) / gammaMom;
  
  std::uniform_real_distribution<double> rndmEngine;
                                           
    constexpr float    iZet        = 13.;
    const double deltaFactor = 136. / pow(Z, oneOverThree) * eps0;
    const double deltaMin    = 4. * deltaFactor;
    
    constexpr double logZ13 = log(iZet) / 3.;
    constexpr double  FZ      =8. * logZ13 + 8. * computeCoulombFactor(); 
    constexpr double  deltaMax = 8.*(logZ13 + computeCoulombFactor());
    
    // compute the limits of eps
    const double epsp        = 0.5 - 0.5*std::sqrt(1. - deltaMin/deltaMax);
    const double epsMin      = std::max(eps0,epsp);
    const double epsRange    = 0.5 - epsMin;
    //
    // sample the energy rate (eps) of the created electron (or positron)
    const double F10 = ScreenFunction1(deltaMin) - FZ;
    const double F20 = ScreenFunction2(deltaMin) - FZ; 
    const double NormF1   = F10 * epsRange * epsRange; 
    const double NormF2   = 1.5 * F20;
    
    // we will need 3 uniform random number for each trial of sampling 
    double greject = 0.;
    double eps;
    do {
      if (NormF1 > rndmEngine(generator) * (NormF1 + NormF2)) {
        eps = 0.5 - epsRange * pow(rndmEngine(generator), oneOverThree);
        const double delta = deltaFactor / (eps * (1. - eps));
          greject = (ScreenFunction1(delta)-FZ)/F10;
      } else {
        eps = epsMin + epsRange * rndmEngine(generator);
        const double delta = deltaFactor / (eps * (1. - eps));
          greject = (ScreenFunction2(delta)-FZ)/F20;
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while (greject < rndmEngine(generator));
    //  end of eps sampling
	return eps;
}

  //~ //
  //~ // select charges randomly
  //~ G4double eTotEnergy, pTotEnergy;
  //~ if (rndmEngine->flat() > 0.5) {
    //~ eTotEnergy = (1.-eps)*gammaMom;
    //~ pTotEnergy = eps*gammaMom;
  //~ } else {
    //~ pTotEnergy = (1.-eps)*gammaMom;
    //~ eTotEnergy = eps*gammaMom;
  //~ }
  //~ //
  //~ // sample pair kinematics
  //~ // 
  //~ const G4double eKinEnergy = std::max(0.,eTotEnergy - CLHEP::electron_mass_c2);
  //~ const G4double pKinEnergy = std::max(0.,pTotEnergy - CLHEP::electron_mass_c2);
  //~ //
  //~ G4ThreeVector eDirection, pDirection;
  //~ //
  //~ GetAngularDistribution()->SamplePairDirections(aDynamicGamma, 
						 //~ eKinEnergy, pKinEnergy, eDirection, pDirection);
  //~ // create G4DynamicParticle object for the particle1
  //~ G4DynamicParticle* aParticle1= new G4DynamicParticle(
                     //~ fTheElectron,eDirection,eKinEnergy);

  //~ // create G4DynamicParticle object for the particle2
  //~ G4DynamicParticle* aParticle2= new G4DynamicParticle(
                     //~ fThePositron,pDirection,pKinEnergy);
  //~ // Fill output vector
  //~ fvect->push_back(aParticle1);
  //~ fvect->push_back(aParticle2);
  //~ // kill incident photon
  //~ fParticleChange->SetProposedKineticEnergy(0.);
  //~ fParticleChange->ProposeTrackStatus(fStopAndKill);
}

  
/// Accoring to the comment in these functions, these are related to F10 and F20
//~ inline G4double G4PairProductionRelModel::ScreenFunction1(G4double ScreenVariable)
//~ // compute the value of the screening function 3*PHI1 - PHI2
//~ {
  //~ return (ScreenVariable > 1.)
    //~ ? 42.24 - 8.368*G4Log(ScreenVariable+0.952)
    //~ : 42.392 - ScreenVariable*(7.796 - 1.961*ScreenVariable);
//~ }

//~ //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//~ inline G4double G4PairProductionRelModel::ScreenFunction2(G4double ScreenVariable)
//~ // compute the value of the screening function 1.5*PHI1 + 0.5*PHI2
//~ {
  //~ return (ScreenVariable > 1.)
    //~ ? 42.24 - 8.368*G4Log(ScreenVariable+0.952)
    //~ : 41.405 - ScreenVariable*(5.828 - 0.8945*ScreenVariable);
//~ }

template <typename generator_t>
Particle::Vector3 
ActsFatras::PhotonConversion::childDirection(generator_t& generator, const ActsFatras::Particle::Vector4& gammaMom4) const
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
void 
ActsFatras::PhotonConversion::recordChilds(generator_t& generator, ActsFatras::Particle& photon,
                                                ActsFatras::Particle::Scalar childEnergy,
                                                const ActsFatras::Particle::Vector3& childDirection,
                                                Acts::PdgParticle pdgProduced) const
{
    // Calculate the child momentum
    const float massChild = findMass(Acts::PdgParticle::eElectron);
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
std::vector<ActsFatras::Particle> 
ActsFatras::PhotonConversion::doConversion(generator_t& generator,
	const ActsFatras::Particle& particle) const {
  const double p = particle.absMomentum();

  // get the energy
  const Particle::Scalar childEnergy = p * childEnergyFraction(generator, p);

  // now get the deflection
  Particle::Vector3 childDir = childDirection(generator, particle.momentum4());
  // verbose output
   return recordChilds(particle,
       childEnergy,
       childDir);
}