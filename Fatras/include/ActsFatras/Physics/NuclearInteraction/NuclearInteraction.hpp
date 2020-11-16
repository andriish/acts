// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Physics/NuclearInteraction/Parameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/Helpers.hpp"
//~ #include "Listener.hpp"

#include <vector>
#include <random>
#include <iterator>
#include <cmath>
#include <any>

namespace ActsFatras {

/// @brief This class provides a parametrised nuclear interaction. The thereby required parametrisation needs to be set and is not provided by default.
struct NuclearInteraction {
	/// The storage of the parameterisation
  detail::MultiParticleParametrisation const* multiParticleParameterisation = nullptr;
  
  //~ Listener listener;
  
  /// @brief Main call operator
  ///
  /// @tparam generator_t The random number generator type
  /// @param [in, out] generator The random number generator
  /// @param [in] slab The material
  /// @param [in, out] particle The ingoing particle
  ///
  /// @return Vector containing the produced secondaries
  template <typename generator_t>
  std::vector<Particle> operator()(generator_t& generator,
                                     const Acts::MaterialSlab& /*slab*/,
                                     Particle& particle) const { return {}; }
	//~ if(multiParticleParameterisation == nullptr)
	//~ {
		//~ return {};
	//~ }
	
	//~ std::uniform_real_distribution<double> uniformDistribution {0., 1.};
	//~ if(particle.pathLimitL0() == std::numeric_limits<Particle::Scalar>::max())
	//~ {
		//~ std::cout << "Need a limit for " << particle.pdg() << std::endl;
		//~ if(multiParticleParameterisation->find(particle.pdg()) != multiParticleParameterisation->end())
		//~ {
			//~ const auto& distribution = multiParticleParameterisation->at(particle.pdg()).at(particle.absMomentum()).nuclearInteractionProbability;
			//~ particle.setMaterialLimits(particle.pathLimitX0(), sampleContinuousValues(uniformDistribution(generator), distribution));
		//~ }
	//~ }
	
	//~ // Test whether enough material was passed for a nuclear interaction
    //~ if(particle.pathInL0() >= particle.pathLimitL0())
    //~ {
		//~ // Fast exit if there is no parametrisation
		//~ if(multiParticleParameterisation->empty())
			//~ return {};
		//~ const auto parametrisationIterator = multiParticleParameterisation->find(particle.pdg());
		//~ if(parametrisationIterator == multiParticleParameterisation->end())
			//~ return {};
		//~ const detail::Parametrisation& parametrisation = parametrisationIterator->second;
		
		//~ // Get the parameters
		//~ const detail::Parameters& parameters = findParameters(uniformDistribution(generator), parametrisation, particle.absMomentum());
//~ const Particle::Vector4 mom4 = particle.momentum4();
		//~ // Dice the interaction type
		//~ std::normal_distribution<double> normalDistribution {0., 1.};
		//~ if(softInteraction(normalDistribution(generator), parameters.softInteractionProbability))
		//~ {
			//~ // Get the final state multiplicity
			//~ const unsigned int multiplicity = finalStateMultiplicity(uniformDistribution(generator), parameters.softMultiplicity); // TODO: test that 0 is forbidden
			//~ // Get the particle content
			//~ const std::vector<int> pdgIds = samplePdgIds(generator, parameters.pdgMap, multiplicity, particle.pdg(), true); // TODO: treat soft interactions
			
			//~ // Get the kinematics
			//~ switch(multiplicity)
			//~ {
				//~ case 1:
				//~ {
					//~ const auto kinematics = sampleKinematics(generator, parameters, 1, true);
					//~ return convertParametersToParticles(generator, pdgIds, kinematics.first, kinematics.second, particle);
				//~ }
				//~ case 2:
				//~ {
					//~ const auto kinematics = sampleKinematics(generator, parameters, 2, true);
					//~ return convertParametersToParticles(generator, pdgIds, kinematics.first, kinematics.second, particle);
				//~ }
				//~ case 3:
				//~ {
					//~ const auto kinematics = sampleKinematics(generator, parameters, 3, true);
					//~ return convertParametersToParticles(generator, pdgIds, kinematics.first, kinematics.second, particle);
				//~ }
				//~ case 4:
				//~ {
					//~ const auto kinematics = sampleKinematics(generator, parameters, 4, true);
					//~ return convertParametersToParticles(generator, pdgIds, kinematics.first, kinematics.second, particle);
				//~ }
				//~ case 5:
				//~ {
					//~ const auto kinematics = sampleKinematics(generator, parameters, 5, true);
					//~ return convertParametersToParticles(generator, pdgIds, kinematics.first, kinematics.second, particle);
				//~ }
			//~ }
			//~ // TODO: Get the initial particle from this set and split the vector accordingly
		//~ } else {
			//~ // Get the final state multiplicity
			//~ const unsigned int multiplicity = finalStateMultiplicity(uniformDistribution(generator), parameters.hardMultiplicity);
			//~ // Get the particle content
			//~ const std::vector<int> pdgIds = samplePdgIds(generator, parameters.pdgMap, multiplicity, particle.pdg(), false);
			
			//~ switch(multiplicity)
			//~ {
				//~ case 0:
				//~ {
					//~ particle.setAbsMomentum(0);
					//~ return {};
				//~ }
				//~ case 1:
				//~ {
					//~ const auto kinematics = sampleKinematics<1>(generator, parameters, false);
					//~ const auto finalStateParticle = convertParametersToParticles<1>(generator, pdgIds, kinematics.first, kinematics.second, particle);
					//~ particle.setAbsMomentum(0);
					//~ return finalStateParticle;
				//~ }
				//~ case 2:
				//~ {
					//~ const auto kinematics = sampleKinematics<2>(generator, parameters, false);
					//~ const auto finalStateParticle = convertParametersToParticles<2>(generator, pdgIds, kinematics.first, kinematics.second, particle);
					//~ particle.setAbsMomentum(0);
					//~ return finalStateParticle;
				//~ }
				//~ case 3:
				//~ {
					//~ const auto kinematics = sampleKinematics<3>(generator, parameters, false);
					//~ const auto finalStateParticle = convertParametersToParticles<3>(generator, pdgIds, kinematics.first, kinematics.second, particle);
					//~ particle.setAbsMomentum(0);
					//~ return finalStateParticle;
				//~ }
				//~ case 4:
				//~ {
					//~ const auto kinematics = sampleKinematics<4>(generator, parameters, false);
					//~ const auto finalStateParticle = convertParametersToParticles<4>(generator, pdgIds, kinematics.first, kinematics.second, particle);
					//~ particle.setAbsMomentum(0);
					//~ return finalStateParticle;
				//~ }
				//~ case 5:
				//~ {
					//~ const auto kinematics = sampleKinematics<5>(generator, parameters, false);
					//~ const auto finalStateParticle = convertParametersToParticles<5>(generator, pdgIds, kinematics.first, kinematics.second, particle);
					//~ particle.setAbsMomentum(0);
					//~ return finalStateParticle;
				//~ }
			//~ }
		//~ }
	//~ }
    //~ // Generates no new particles
    //~ return {};
  //~ }
  
  private:
    /// @brief Retrieves the parametrisation for the particle
    ///
    /// @param [in] rnd Random number
    /// @param [in] particle The current particle
    ///
    /// @return The parametrisation
    const detail::Parameters& findParameters(double rnd, const detail::Parametrisation& parametrisation, float particleMomentum) const;
    
    /// @brief Estimates the interaction type
    ///
    /// @param [in] rnd Random number
    /// @param [in] probability The probability for a soft interaction
    ///
    /// @return True if a soft interaction occurs
    bool softInteraction(double rnd, float probability) const { return rnd <= probability; }
    
    /// @brief Evaluates the multiplicity of the final state
    ///
    /// @param [in] rnd Random number
    /// @param [in] distribution The multiplicity distribution
    ///
    /// @return The final state multiplicity
    unsigned int finalStateMultiplicity(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const;
    
    /// @brief Evaluates the final state PDG IDs
    ///
    /// @tparam generator_t The random number generator type
    /// @param [in, out] generator The random number generator
    /// @param [in] multiplicity The final state multiplicity
    ///
    /// @return Vector containing the PDG IDs
    template <typename generator_t>
    std::vector<int> samplePdgIds(generator_t& generator, const detail::Parameters::PdgMap& pdgMap, unsigned int multiplicity, int particlePdg, bool soft) const;

    /// @brief Evaluates the final state invariant masses
    ///
    /// @tparam generator_t The random number generator type
    /// @param [in, out] generator The random number generator
    /// @param [in] kinematicParameters Parametrisation of kinematic properties
    ///
    /// @return Vector containing the invariant masses
    template <typename generator_t>
    Acts::ActsVectorXf sampleInvariantMasses(generator_t& generator, const detail::Parameters::ParametersWithFixedMultiplicity& parametrisation) const;
      
    /// @brief Evaluates the final state momenta
    ///
    /// @tparam generator_t The random number generator type
    /// @param [in, out] The random number generator
    /// @param [in] kinematicParameters Parametrisation of kinematic properties
    ///
    /// @return Vector containing the momenta
    template <typename generator_t>
    Acts::ActsVectorXf sampleMomenta(generator_t& generator, const detail::Parameters::ParametersWithFixedMultiplicity& parametrisation,
		float initialMomentum) const;
	
	template<typename generator_t>
	std::pair<Acts::ActsVectorXf, Acts::ActsVectorXf> sampleKinematics(generator_t& generator, 
					const detail::Parameters& parameters, unsigned int multiplicity, bool soft) const
	{
		// TODO: let the any cast disappear
		auto const* kinematicParameters = std::any_cast<const detail::Parameters::ParametersWithFixedMultiplicity>
			(soft ? (&parameters.softKinematicParameters[multiplicity]) : (&parameters.hardKinematicParameters[multiplicity]));
		const auto invariantMasses = sampleInvariantMasses(generator, *kinematicParameters);
		auto momenta = sampleMomenta(generator, *kinematicParameters, parameters.momentum);
		while(!match(momenta, invariantMasses, parameters.momentum)) {
			momenta = sampleMomenta(generator, *kinematicParameters, parameters.momentum);
		}
		return std::make_pair(momenta, invariantMasses);
	}
    
bool
match(const Acts::ActsVectorXf& momenta, const Acts::ActsVectorXf& invariantMasses, float initialMomentum) const;

std::pair<ActsFatras::Particle::Scalar, ActsFatras::Particle::Scalar>
globalAngle(ActsFatras::Particle::Scalar phi1, ActsFatras::Particle::Scalar theta1, float phi2, float theta2) const;
  
    float pathLimitL0(double rnd, detail::Parameters::CumulativeDistribution& distribution) const;
    
    /// @brief Converter from sampled numbers to a vector of particles
    ///
    /// @param [in] pdgId The PDG IDs
    /// @param [in] momentum The momenta
    /// @param [in] invariantMass The invariant masses
    ///
    /// @return Vector containing the final state particles
    template <typename generator_t>
	std::vector<Particle> convertParametersToParticles(generator_t& generator, const std::vector<int>& pdgId, const Acts::ActsVectorXf& momenta, 
							const Acts::ActsVectorXf& invariantMasses, const Particle& initialParticle) const;
    
    ///Function gets random number rnd in the range [0,1) as argument 
    ///and returns function value according to a histogram distribution
    unsigned int sampleDiscreteValues(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const; // TODO: this should be only used for multiplicity
	
	
	double sampleContinuousValues(double rnd, const detail::Parameters::CumulativeDistribution& distribution, bool interpolate = false) const;
};

	template <typename generator_t>
	Acts::ActsVectorXf
	NuclearInteraction::sampleInvariantMasses(generator_t& generator, const detail::Parameters::ParametersWithFixedMultiplicity& parametrisation) const {
		
	  Acts::ActsVectorXf parameters;
	  const unsigned int size = parametrisation.eigenvaluesInvariantMass.size();
	  
	  for(unsigned int i = 0; i < size; i++) {
		float variance = parametrisation.eigenvaluesInvariantMass[i];
		std::normal_distribution<float> dist{parametrisation.meanInvariantMass[i], sqrtf(variance)};
		parameters[i] = dist(generator);
	  }
	  parameters = parametrisation.eigenvectorsInvariantMass * parameters;
	  
		for(unsigned int i = 0; i < size; i++)
		{
			const double cdf = (std::erff(parameters[i]) + 1) * 0.5;
			parameters[i] = sampleContinuousValues(cdf, parametrisation.invariantMassDistributions[i]);
		}
		return parameters;
	}
	
	template <typename generator_t>
	Acts::ActsVectorXf
	NuclearInteraction::sampleMomenta(generator_t& generator, const detail::Parameters::ParametersWithFixedMultiplicity& parametrisation, float initialMomentum) const {
		
	  Acts::ActsVectorXf parameters;
	  const unsigned int size = parametrisation.eigenvaluesMomentum.size();
	
	  for(unsigned int i = 0; i < size + 1; i++) {
		float variance = parametrisation.eigenvaluesMomentum[i];
		std::normal_distribution<float> dist{parametrisation.meanMomentum[i], sqrtf(variance)};
		parameters[i] = dist(generator);
	  }
	  parameters = parametrisation.eigenvectorsMomentum * parameters;
	  
		for(unsigned int i = 0; i < size + 1; i++)
		{
			const float cdf = (std::erff(parameters[i]) + 1) * 0.5;
			parameters[i] = sampleContinuousValues(cdf, parametrisation.momentumDistributions[i]);
		}
		
		Acts::ActsVectorXf momenta = parameters.template head(size);
		const float sum = momenta.sum();
		const float scale = parameters.template tail<1>()(0,0) / sum;
		momenta *= scale * initialMomentum;
		
		return momenta;
	}
	
	
template <typename generator_t>
std::vector<int> NuclearInteraction::samplePdgIds(generator_t& generator, const detail::Parameters::PdgMap& pdgMap, 
	unsigned int multiplicity, int particlePdg, bool soft) const {
	
	if(multiplicity == 0)
		return {};
	
	std::vector<int> pdgIds;
	pdgIds.reserve(multiplicity);
	
	std::uniform_real_distribution<float> uniformDistribution {0., 1.};
	
	auto citProducer = pdgMap.cbegin();
	while(citProducer->first != particlePdg && citProducer != pdgMap.end())
		citProducer++;
		
	if(citProducer == pdgMap.end()) // TODO: if a parametrisation is available then this should never occur
	{
		return {}; // TODO: handle this error
	}
	const std::vector<std::pair<int, float>>& mapInitial = citProducer->second;
	
	if(soft)
		// Store the initial particle if the interaction is soft
		pdgIds.push_back(particlePdg);
	else
	{
		// Otherwise dice the particle
		const float rndInitial = uniformDistribution(generator);
		pdgIds.push_back(std::lower_bound(mapInitial.begin(), mapInitial.end(), rndInitial, [](const std::pair<int, float>& element, float random)
			{ return element.second < random; })->second);
	}
	
	for(unsigned int i = 1; i < multiplicity; i++)
	{ // TODO: could there be the case, that a map doesn't exist / is empty? - maybe set all keys to avoid this
		citProducer = pdgMap.cbegin();
		while(citProducer->first != pdgIds[i - 1] && citProducer != pdgMap.end())
			citProducer++;
		
		const std::vector<std::pair<int, float>>& map = citProducer->second;
		const float rnd = uniformDistribution(generator);
		pdgIds.push_back(
				std::lower_bound(map.begin(), map.end(), rnd, [](const std::pair<int, float>& element, float random){ return element.second < random; })->second);
	}
	return {};
}

template <typename generator_t>
std::vector<Particle> NuclearInteraction::convertParametersToParticles(generator_t& generator, const std::vector<int>& pdgId, const Acts::ActsVectorXf& momenta, 
							const Acts::ActsVectorXf& invariantMasses, const Particle& initialParticle) const {
  std::uniform_real_distribution<double> uniformDistribution {0., 1.};
  std::vector<Particle> result;	
// TODO: handle soft case
// TODO: dice soft particle L0 again
  const auto initialMomentum = initialParticle.absMomentum();
  const auto& initialDirection = initialParticle.unitDirection();
  const double phi = Acts::VectorHelpers::phi(initialDirection);
  const double theta = Acts::VectorHelpers::theta(initialDirection);
  const unsigned int size = momenta.size();
  for(unsigned int i = 0; i < size; i++)
  {
		const float momentum = momenta[i];
		const float invariantMass = invariantMasses[i];	
		const float p1p2 = 2. * momentum * initialMomentum;
		const float costheta = 1. - invariantMass * invariantMass / p1p2;
		
	  const auto phiTheta = globalAngle(phi, theta, uniformDistribution(generator) * 2. * M_PI, std::acos(costheta));
	  const auto direction = Acts::makeDirectionUnitFromPhiTheta(phiTheta.first, phiTheta.second);
	  Particle p = Particle(Barcode(), static_cast<Acts::PdgParticle>(pdgId[i]));
	  p.setProcess(ProcessType::eNuclearInteraction).setPosition4(initialParticle.position4())
		.setAbsMomentum(momentum).setDirection(direction);
	
	auto cit = multiParticleParameterisation->begin();
	while(cit->first != p.pdg() && cit != multiParticleParameterisation->end())
		cit++;
		
	if(cit != multiParticleParameterisation->end())
	{
		const auto& distribution = findParameters(uniformDistribution(generator), cit->second, p.absMomentum()).nuclearInteractionProbability;
		p.setMaterialLimits(p.pathLimitX0(), sampleContinuousValues(uniformDistribution(generator), distribution));
	}
  }
  					
  return result;
}
}  // namespace ActsFatras
