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

#include <vector>
#include <random>
#include <iterator>
#include <cmath>

namespace ActsFatras {

/// @brief This class provides a parametrised nuclear interaction. The thereby required parametrisation needs to be set and is not provided by default.
struct NuclearInteraction {
	/// The storage of the parameterisation
  detail::MultiParticleParametrisation multiParticleParameterisation; // TODO: this should depend on the particle
  
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
                                     const Acts::MaterialSlab& slab,
                                     Particle& particle) const {
	// Test whether enough material was passed for a nuclear interaction
    if(particle.pathInL0() >= particle.pathLimitL0())
    {
		// Fast exit if there is no parametrisation
		if(multiParticleParameterisation.empty())
			return {};
		const auto parametrisationIterator = multiParticleParameterisation.find(particle.pdg());
		if(parametrisationIterator == multiParticleParameterisation.end())
			return {};
		const detail::Parametrisation& parametrisation = *parametrisationIterator;
		
		std::uniform_real_distribution<double> uniformDistribution {0., 1.};
		std::normal_distribution<double> normalDistribution {0., 1.};
		
		// Get the parameters
		const detail::Parameters& parameters = findParameters(uniformDistribution(generator), particle);

		// Dice the interaction type
		if(softInteraction(normalDistribution(generator), parameters.softInteractionProbability))
		{
			// Get the final state multiplicity
			const unsigned int multiplicity = finalStateMultiplicity(uniformDistribution(generator), parameters.softMultiplicity); // TODO: test that 0 is forbidden
			// Get the particle content
			const std::vector<int> pdgIds = samplePdgIds(generator, parameters.pdgMap, multiplicity, particle.pdg()); // TODO: treat soft interactions
			
			// Get the kinematics
			const auto& kinematicParameters = parameters.softKinematicParameters[multiplicity];
			const auto invariantMasses = sampleInvariantMasses(generator, kinematicParameters);
			auto momenta = sampleMomenta(generator, kinematicParameters, parameters.momentum);
			while(!match(momenta, invariantMasses, parameters.momentum)) {
				momenta = sampleMomenta(generator, kinematicParameters, parameters.momentum);
			}
			
			// Build and return particles
			return convertParametersToParticles(pdgIds, momenta, invariantMasses, particle); // TODO: Get the initial particle from this set and split the vector accordingly
		} else {
			// Get the final state multiplicity
			const unsigned int multiplicity = finalStateMultiplicity(uniformDistribution(generator), parameters.hardMultiplicity);
			// Get the particle content
			const std::vector<int> pdgIds = samplePdgIds(generator, parameters.pdgMap, multiplicity, particle.pdg());
			
			// Get the kinematics
			const auto& kinematicParameters = parameters.hardKinematicParameters[multiplicity];
			const auto invariantMasses = sampleInvariantMasses(generator, kinematicParameters);
			auto momenta = sampleMomenta(generator, kinematicParameters, parameters.momentum);
			while(!match(momenta, invariantMasses, parameters.momentum)) {
				momenta = sampleMomenta(generator, kinematicParameters, parameters.momentum);
			}
			
			// Build and return particles
			return convertParametersToParticles(pdgIds, momenta, invariantMasses, particle);
		}
	}
    // Generates no new particles
    return {};
  }
  
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
    std::vector<int> samplePdgIds(generator_t& generator, const detail::Parameters::PdgMap& pdgMap, unsigned int multiplicity, int particlePdg) const;

    /// @brief Evaluates the final state invariant masses
    ///
    /// @tparam generator_t The random number generator type
    /// @param [in, out] generator The random number generator
    /// @param [in] kinematicParameters Parametrisation of kinematic properties
    ///
    /// @return Vector containing the invariant masses
    template <unsigned int size, typename generator_t>
    Acts::ActsVectorF<size> sampleInvariantMasses(generator_t& generator, const detail::Parameters::ParametersWithFixedMultiplicity<size>& parametrisation) const;
      
    /// @brief Evaluates the final state momenta
    ///
    /// @tparam generator_t The random number generator type
    /// @param [in, out] The random number generator
    /// @param [in] kinematicParameters Parametrisation of kinematic properties
    ///
    /// @return Vector containing the momenta
    template <unsigned int size, typename generator_t>
    Acts::ActsVectorF<size - 1> sampleMomenta(generator_t& generator, const detail::Parameters::ParametersWithFixedMultiplicity<size>& parametrisation,
		float initialMomentum) const;
    
template<unsigned int size>
bool
match(const Acts::ActsVectorF<size>& momenta, const Acts::ActsVectorF<size>& invariantMasses, float initialMomentum) const;

std::pair<float, float>
globalAngle(float phi1, float theta1, float phi2, float theta2) const;
  
    float pathLimitL0(double rnd, detail::Parameters::CumulativeDistribution& distribution) const;
    
    /// @brief Converter from sampled numbers to a vector of particles
    ///
    /// @param [in] pdgId The PDG IDs
    /// @param [in] momentum The momenta
    /// @param [in] invariantMass The invariant masses
    ///
    /// @return Vector containing the final state particles
    template <unsigned int size>
	std::vector<Particle> convertParametersToParticles(const std::vector<int>& pdgId, const Acts::ActsVectorF<size>& momentum, 
							const Acts::ActsVectorF<size>& invariantMass, const Particle& initialParticle) const;
// TODO: the conversion should set the limit_l0	
    ///Function gets random number rnd in the range [0,1) as argument 
    ///and returns function value according to a histogram distribution
    unsigned int sampleDiscreteValues(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const; // TODO: this should be only used for multiplicity
	
	
	double sampleContinuousValues(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const;
};

	template <unsigned int size, typename generator_t>
	Acts::ActsVectorF<size>
	NuclearInteraction::sampleInvariantMasses(generator_t& generator, const detail::Parameters::ParametersWithFixedMultiplicity<size>& parametrisation) const {
		
	  Acts::ActsVectorF<size> parameters;
		
	  for(unsigned int i = 0; i < size; i++) {
		float variance = parametrisation.eigenvaluesInvariantMass[i];
		std::normal_distribution<float> dist{parametrisation.meanInvariantMass[i], sqrt(variance)};
		parameters[i] = dist(generator);
	  }
	  parameters = parametrisation.eigenvectorsInvariantMass * parameters;
	  
		for(int i = 0; i < size; i++)
		{
			const double cdf = (std::erff(parameters[i]) + 1) * 0.5;
			parameters[i] = sampleContinuousValues(cdf, parametrisation.invariantMassDistributions[i]);
		}
		return parameters;
	}
	
	template <unsigned int size, typename generator_t>
	Acts::ActsVectorF<size - 1>
	NuclearInteraction::sampleMomenta(generator_t& generator, const detail::Parameters::ParametersWithFixedMultiplicity<size>& parametrisation, float initialMomentum) const {
		
	  Acts::ActsVectorF<size> parameters;
		
	  for(unsigned int i = 0; i < size; i++) {
		float variance = parametrisation.eigenvaluesMomentum[i];
		std::normal_distribution<float> dist{parametrisation.meanMomentum[i], sqrt(variance)};
		parameters[i] = dist(generator);
	  }
	  parameters = parametrisation.eigenvectorsMomentum * parameters;
	  
		for(int i = 0; i < size; i++)
		{
			const float cdf = (std::erff(parameters[i]) + 1) * 0.5;
			parameters[i] = sampleContinuousValues(cdf, parametrisation.momentumDistributions[i]);
		}
		
		Acts::ActsVectorF<size - 1> momenta = parameters.template head<size - 1>();
		const float sum = momenta.sum();
		const float scale = parameters.template tail<1>() / sum;
		momenta *= scale * initialMomentum;
		
		return momenta;
	}
	
	
template <typename generator_t>
std::vector<int> NuclearInteraction::samplePdgIds(generator_t& generator, const detail::Parameters::PdgMap& pdgMap, 
	unsigned int multiplicity, int particlePdg) const {
	std::vector<int> pdgIds; // TODO: handle soft interactions
	pdgIds.reserve(multiplicity);
	
	std::uniform_real_distribution<float> uniformDistribution {0., 1.};
	
	if(pdgMap.find(particlePdg) == pdgMap.end()) // TODO: if a parametrisation is available then this should never occur
	{
		return {}; // TODO: handle this error
	}
	const std::unordered_map<int, float>& mapInitial = pdgMap.at(particlePdg);
	const float rndInitial = uniformDistribition(generator);
	pdgIds.push_back(std::lower_bound(mapInitial.begin(), mapInitial.end(), rndInitial, [](const std::pair<int, float>& element, float random)
		{ return element.second < random; }));
	
	for(unsigned int i = 1; i < multiplicity; i++)
	{ // TODO: could there be the case, that a map doesn't exist / is empty? - maybe set all keys to avoid this
		const std::unordered_map<int, float>& map = pdgMap.at(pdgIds[i - 1]);
		const float rnd = uniformDistribition(generator);
		pdgIds.push_back(std::lower_bound(map.begin(), map.end(), rnd, [](const std::pair<int, float>& element, float random){ return element.second < random; }));
	}
	return {};
}

template<unsigned int size>
bool
NuclearInteraction::match(const Acts::ActsVectorF<size>& momenta, const Acts::ActsVectorF<size>& invariantMasses, float initialMomentum) const
{
	for(unsigned int i = 0; i < size; i++)
	{
		const float momentum = momenta[i];
		const float invariantMass = invariantMasses[i];
		
		const double p1p2 = 2. * momentum * initialMomentum;
		const double costheta = 1. - invariantMass * invariantMass / p1p2;

		if(std::abs(costheta) > 1)
		{
			return false;
		}
	}
	return true;
}

template <unsigned int size>
std::vector<Particle> NuclearInteraction::convertParametersToParticles(const std::vector<int>& pdgId, const Acts::ActsVectorF<size>& momentum, 
							const Acts::ActsVectorF<size>& invariantMass, const Particle& initialParticle) const {
  std::vector<Particle> result;	

  for(unsigned int i = 0; i < size; i++)
  {
	  Particle p = Particle(Barcode(), pdgId[i]).setProcess(ProcessType::eNuclearInteraction).setPosition4(initialParticle.position4()).absMomentum(momentum[i]);
  }

  //~ Vector3 m_unitDirection = Vector3::UnitZ();


  //~ Scalar m_limitX0 = std::numeric_limits<Scalar>::max();
  //~ Scalar m_limitL0 = std::numeric_limits<Scalar>::max();
  					
  return result;
}
}  // namespace ActsFatras
