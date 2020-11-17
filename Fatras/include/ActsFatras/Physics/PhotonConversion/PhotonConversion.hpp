// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Gaudi & Athena
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "AthenaKernel/IAtRndmGenSvc.h"
// Trk
#include "TrkEventPrimitives/PropDirection.h"
#include "TrkEventPrimitives/ParticleHypothesis.h"
#include "TrkExUtils/MaterialUpdateMode.h"
#include "TrkDetDescrUtils/GeometrySignature.h" 
// ISF
#include "ISF_Event/ITruthIncident.h"
// Fatras
#include "ISF_FatrasInterfaces/IPhotonConversionTool.h"
// Barcode
#include "BarcodeEvent/PhysicsProcessCode.h"



#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"

namespace Trk {
  class Layer;
  class CylinderVolumeBounds;
  class PdgToParticleHypothesis;
  class TrackingGeometry;
  class ITrackingGeometrySvc;
}

namespace ActsFatras {

  
  /** @class PhotonConversionTool
      
     The photon conversion tool, to be called by the MaterialUpdator.
  
     @author Sarka.Todorova@cern.ch
  */
   
  class PhotonConversionTool {
  public:          
      /** interface for processing of the pair production */
      bool pairProduction(const Trk::MaterialProperties& mprop,
			  double pathCorrection,
			  double p) const;
      
      /** interface for processing of the presampled pair production */      
      /** interface for processing of the presampled nuclear interactions on layer*/
      std::vector<Particle> doConversion(double time, const Trk::NeutralParameters& parm) const;


    
   private:
      /** record childs - create interface already for mu-/mu+ pair production*/
      void recordChilds(const Particle& photon, 
			double childEnergy,
			const Acts::Vector3D& childDirection,
			Acts::PdgParticle = 11) const;
	
      /** simulate the child energy */
      double childEnergyFraction(double gammaMom) const;
      
      /** simulate the one child direction - the second one is given clear then*/
      Amg::Vector3D childDirection(const Amg::Vector3D& gammaMom,
					  double childE) const;
      
      /** helper functions for the Phi1/phi2 */
      double phi1(double delta) const;
      
      /** helper functions for the Phi1/phi2 */
      double phi2(double delta) const;
      
      /** MCTruth process code for TruthIncidents created by this tool */
      Barcode::PhysicsProcessCode                  m_processCode = 14; // TODO: Will become part of the process enum
      
      /** The cut from which on the child products are followed */
      double                                       m_minChildEnergy = 50. * Acts::UnitConstants::MeV;
      double                                       m_childEnergyScaleFactor = 2.;
      double                                       m_conversionProbScaleFactor = 0.98;
            
      /** struct of Particle Masses */
      static Trk::PdgToParticleHypothesis          s_pdgToHypo;
      /** Inverse fine structure constant */
      constexpr double                                s_alphaEM = 1. / 137.;
      constexpr double                                s_oneOverThree = 1. / 3.;
   };
  

inline double PhotonConversionTool::phi1(double delta) const {
  if (delta <= 1.)
     return 20.867 - 3.242 * delta  + 0.625*delta*delta;
  else
    return 21.12 - 4.184*log(delta+0.952);
}

inline double PhotonConversionTool::phi2(double delta) const {
  if (delta <= 1.)
     return 20.209 - 1.930 * delta  - 0.086*delta*delta;
   return 21.12 - 4.184*log(delta+0.952);
}

}



///////////////////////////////////////////////////////////////////
// PhotonConversionTool.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// class header
#include "PhotonConversionTool.h"

// Gaudi Kernel
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
// ISF includes
#include "ISF_Event/ISFParticle.h"
#include "ISF_Event/ISFParticleVector.h"
#include "ISF_Event/ParticleClipboard.h"
#include "ISF_Event/ParticleUserInformation.h"
#include "ISF_Event/ISFTruthIncident.h"
// Trk inlcude
#include "TrkEventPrimitives/PdgToParticleHypothesis.h"
#include "TrkExInterfaces/IEnergyLossUpdator.h"
#include "TrkExInterfaces/ITimedExtrapolator.h"
#include "TrkExInterfaces/IMultipleScatteringUpdator.h"
#include "TrkSurfaces/Surface.h"
#include "TrkGeometry/Layer.h"
#include "TrkGeometry/MaterialProperties.h"
#include "TrkVolumes/CylinderVolumeBounds.h"
#include "TrkNeutralParameters/NeutralParameters.h"
// CLHEP
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzVector.h"
// STD
#include <math.h>

// temporary
#include "TrkDetDescrInterfaces/ITrackingGeometrySvc.h"
#include "TrkGeometry/TrackingGeometry.h"
#include "TrkGeometry/TrackingVolume.h"

// statics doubles 
Trk::PdgToParticleHypothesis  ActsFatras::PhotonConversionTool::s_pdgToHypo;


// ----------------- private helper methods ----------------------------------------------------
template <typename generator_t>
void ActsFatras::PhotonConversionTool::recordChilds(generator_t& generator, ActsFatras::Particle& photon,
                                                double childEnergy,
                                                const Acts::Vector3D& childDirection,
                                                Acts::PdgParticle pdgProduced) const
{
    // Calculate the child momentum
    const double massChild = findMass(pdgProduced);
    const double momentum1 = sqrt(childEnergy * childEnergy - massChild * massChild);    

    // now properly : energy-momentum conservation
    const Particle::Vector3 vtmp = photom.momentum4().template segment<3>(Acts::eDir0) - momentum1 * childDirection;
    const double momentum2 = vtmp.norm();

    // charge sampling
    std::uniform_real_distribution<double> uniformDistribution {0., 1.};
    Particle::Scalar charge1;
    Particle::Scalar charge2;
    charge1 = charge2 = 0.;
    if (CLHEP::RandFlat::shoot(m_randomEngine)>0.5) {
      charge1 = -1.;
      charge2 =  1.;
    }
    else {
      charge1 =  1.;
      charge2 = -1.;
    }

    // add the new secondary states to the ISF particle stack
    int    pdg1  = s_pdgToHypo.convert(childType, charge1, false);
    int    pdg2  = s_pdgToHypo.convert(childType, charge2, false);

    // remove soft children
    int nchild = 0;
    if ( p1 > m_minChildEnergy ) nchild++;
    if ( p2 > m_minChildEnergy ) nchild++;

    ISF::ISFParticleVector children(nchild);

    int ichild = 0;
    if (  p1 > m_minChildEnergy ) {
      ISF::ISFParticle* ch1 = new ISF::ISFParticle( vertex,
                                               p1*childDirection,
                                               mass,
                                               charge1,
                                               pdg1,
                                               time,
                                               *parent );
      children[ichild] = ch1;
      ichild++;
    }

    if (  p2 > m_minChildEnergy ) {
      ISF::ISFParticle* ch2  = new ISF::ISFParticle( vertex,
                                               p2*childDirection,
                                               mass,
                                               charge2,
                                               pdg2,
                                               time,
                                               *parent );
      children[ichild] = ch2;
    }

    // register TruthIncident
    ISF::ISFTruthIncident truth( const_cast<ISF::ISFParticle&>(*parent),
                                 children,
                                 m_processCode,
                                 parent->nextGeoID(),
                                 ISF::fKillsPrimary );
    m_truthRecordSvc->registerTruthIncident( truth);

}

bool ActsFatras::PhotonConversionTool::pairProduction(const Trk::MaterialProperties& mprop,
                                                 double pathCorrection,
                                                 double p) const
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

 double  p0  =       -7.01612e-03;
 double  p1  =        7.69040e-02;
 double  p2  =       -6.07682e-01;
 // calculate xi
 double xi = p0 + p1*pow(p/1000.,p2);
 // now calculate what's left
 double attenuation = exp( -7.777e-01*pathCorrection*mprop.thicknessInX0()*(1.-xi) );

  return (m_conversionProbScaleFactor*CLHEP::RandFlat::shoot(m_randomEngine) > attenuation) ? true : false;
  // TODO: Transform this probability to a sample in X_0

}


double ActsFatras::PhotonConversionTool::childEnergyFraction(double gammaMom) const {

  // the fraction
  double epsilon0      = findMass(Trk::electron)/gammaMom;
  // some needed manipolations
  double Z             = 13.; //mprop.averageZ(); // TODO: the other part was never used
  double oneOverZpow   = 1./pow(Z,s_oneOverThree);
  double alphaZsquare  = (s_alpha*s_alpha*Z*Z);
  // now f(Z) - from Z and s_alpha
  double fZ            = alphaZsquare*(1./(1.+alphaZsquare)+0.20206-0.0369*alphaZsquare+0.0083*alphaZsquare*alphaZsquare);
  // delta_max
  double deltaMax      = exp((42.24-fZ)*.1195)-0.952;
  // delta_min
  double deltaMin      = 4.*epsilon0*136.*oneOverZpow; 
  // the minimum fraction
  double epsilon1      = 0.5-0.5*sqrt(1.-deltaMin/deltaMax);
  double epsilonMin    = epsilon1 > epsilon0 ? epsilon1 : epsilon0;
  // calculate phi1 / phi2 - calculate from deltaMin
  double Phi1          = phi1(deltaMin);
  double Phi2          = phi2(deltaMin);
  // then calculate F10/F20
  double F10           = 3.*Phi1 - Phi2 - fZ;
  double F20           = 1.5*Phi1 - 0.5*Phi2 - fZ;
  // and finally calucate N1, N2
  double N1            = (0.25-epsilonMin+epsilonMin*epsilonMin)*F10;
  double N2            = 1.5*F20;
  // ------------ decide wich one to take 
  if ( N1/(N1+N2) < CLHEP::RandFlat::shoot(m_randomEngine) ) {  
    // sample from f1,g1 distribution
    for ( ; ; ){
      double epsilon = 0.5 - (0.5 - epsilonMin)*pow(CLHEP::RandFlat::shoot(m_randomEngine),s_oneOverThree);
      // prepare for the rejection check
      double delta   = 136.*epsilon0*oneOverZpow/(epsilon-epsilon*epsilon);
      double F1 = 3.*phi1(delta)-phi2(delta)-fZ;   
      // reject ? - or redo the exercise 
      if (F1/F10 > CLHEP::RandFlat::shoot(m_randomEngine)) return m_childEnergyScaleFactor*epsilon;
    }
  } else {
    // sample from f2,g2 distribution
    for ( ; ; ){
      double epsilon = epsilonMin + (0.5-epsilonMin)*CLHEP::RandFlat::shoot(m_randomEngine);
      // prepare for the rejection check
      double delta   = 136.*epsilon0*oneOverZpow/(epsilon-epsilon*epsilon);
      double F2 = 1.5*phi1(delta)-0.5*phi2(delta)-fZ;   
     // reject ? - or redo the exercise 
     if (F2/F20 > CLHEP::RandFlat::shoot(m_randomEngine)) return m_childEnergyScaleFactor*epsilon;  
    }
  }

}

Amg::Vector3D ActsFatras::PhotonConversionTool::childDirection(const Amg::Vector3D& gammaMom,
                                                                 double childE) const
{
    // --------------------------------------------------
    // Following the Geant4 approximation from L. Urban
    // the azimutal angle
    double psi    =  2.*M_PI*CLHEP::RandFlat::shoot(m_randomEngine);
    
    // the start of the equation
    double theta = findMass(Trk::electron)/childE;
    // follow 
    double a = 0.625; // 5/8
    //double d = 27.;

    double r1 = CLHEP::RandFlat::shoot(m_randomEngine);
    double r2 = CLHEP::RandFlat::shoot(m_randomEngine);
    double r3 = CLHEP::RandFlat::shoot(m_randomEngine);

    double u =  -log(r2*r3)/a;
    
    theta *= (r1 < 0.25 ) ? u : u*s_oneOverThree; // 9./(9.+27) = 0.25

     ATH_MSG_VERBOSE( "[ conv ] Simulated angle to photon    = " << theta << "." );

    // more complex but "more true"
    CLHEP::Hep3Vector gammaMomHep( gammaMom.x(), gammaMom.y(), gammaMom.z() );
    CLHEP::Hep3Vector newDirectionHep(gammaMomHep.unit());
    double x = -newDirectionHep.y();
    double y = newDirectionHep.x();
    double z = 0.;
    // if it runs along the z axis - no good ==> take the x axis
    if (newDirectionHep.z()*newDirectionHep.z() > 0.999999)       
        x = 1.;
    // deflector direction
    CLHEP::Hep3Vector deflectorHep(x,y,z);
    // rotate the new direction for scattering
    newDirectionHep.rotate(theta, deflectorHep);
    // and arbitrarily in psi             
    newDirectionHep.rotate(psi,gammaMomHep);

    // assign the new values
    Amg::Vector3D newDirection( newDirectionHep.x(), newDirectionHep.y(), newDirectionHep.z() );
    return newDirection;

}

/** interface for processing of the presampled nuclear interactions on layer*/
std::vector<ActsFatras::Particle> ActsFatras::PhotonConversionTool::doConversion(
	double time, const Trk::NeutralParameters& parm) const {
  double p = parm.momentum().mag();

  // get the energy
  double childEnergy = p*childEnergyFraction(p);

  // now get the deflection
  Amg::Vector3D childDir(childDirection(parm.momentum(), childEnergy));
  // verbose output
  ATH_MSG_VERBOSE(  "[ conv ] Child energy simulated as : " << childEnergy << " MeV" );
  
	    recordChilds(time,
	       parm.position(),
               parm.momentum().unit(),
	       childEnergy, p,
	       childDir,
	       Trk::electron);
}