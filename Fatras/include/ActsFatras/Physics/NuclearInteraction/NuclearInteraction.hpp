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

#include <vector>
#include <random>

namespace ActsFatras {

/// @brief This class provides a parametrised nuclear interaction. The thereby required parametrisation needs to be set and is not provided by default.
struct NuclearInteraction {
	/// The storage of the parameterisation
  detail::Parametrisation parameterisation; // TODO: this should depend on the particle
  
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
	// Fast exit if there is no parametrisation
	if(parameterisation.empty())
		return {};
	
	std::uniform_real_distribution<double> uniformDistribution {0., 1.};
    std::normal_distribution<double> normalDistribution {0., 1.};
    
	// Get the parameters
	const detail::Parameters& parameters = findParameters(uniformDistribution(generator), particle);

	// Dice whether there is a nuclear interaction
    if(nuclearInteraction(normalDistribution(generator), parameters))
    {
		// Dice the interaction type
		if(softInteraction(normalDistribution(generator), parameters.softInteractionProbability))
		{
			// Get the final state multiplicity
			const unsigned int multiplicity = finalStateMultiplicity(uniformDistribution(generator), parameters.softMultiplicity);
			// Get the particle content
			const std::vector<int> pdgIds = samplePdgIds(generator, parameters.pdgMap, multiplicity);
			
			// Get the kinematics
			const auto& kinematicParameters = parameters.softKinematicParameters[multiplicity];
			const std::vector<float> momenta = sampleMomenta(generator, kinematicParmaeters);
			const std::vector<float> invariantMasses = sampleInvariantMasses(generator, kinematicParmaeters);
			
			// Build and return particles
			return convertParametersToParticles(pdgIds, momenta, invariantMasses); // TODO: Get the initial particle from this set and split the vector accordingly
		} else {
			// Get the final state multiplicity
			const unsigned int multiplicity = finalStateMultiplicity(uniformDistribution(generator), parameters.hardMultiplicity);
			// Get the particle content
			const std::vector<int> pdgIds = samplePdgIds(generator, parameters.pdgMap, multiplicity);
			
			// Get the kinematics
			const auto& kinematicParameters = parameters.hardKinematicParameters[multiplicity];
			const std::vector<float> momenta = sampleMomenta(generator, kinematicParmaeters);
			const std::vector<float> invariantMasses = sampleInvariantMasses(generator, kinematicParmaeters);
			
			// Build and return particles
			return convertParametersToParticles(pdgIds, momenta, invariantMasses);
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
    const detail::Parameters& findParameters(double rnd, const Particle& particle) const;
    
    /// @brief Estimates whether a nuclear interaction occurs
    ///
    /// @param [in] rnd Random number
    /// @param [in] distribution The nuclear interaction probability distribution
    ///
    /// @return True if a nuclear interaction occurs
    bool nuclearInteraction(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const;
    
    /// @brief Estimates the interaction type
    ///
    /// @param [in] rnd Random number
    /// @param [in] probability The probability for a soft interaction
    ///
    /// @return True if a soft interaction occurs
    bool softInteraction(double rnd, float probability) const;
    
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
    std::vector<int> samplePdgIds(generator_t& generator, unsigned int multiplicity) const;
    
    /// @brief Evaluates the final state momenta
    ///
    /// @tparam generator_t The random number generator type
    /// @param [in, out] The random number generator
    /// @param [in] kinematicParameters Parametrisation of kinematic properties
    ///
    /// @return Vector containing the momenta
    template <typename generator_t>
    std::vector<float> sampleMomenta(generator_t& generator, const std::any& kinematicParmaeters) const;
    
    /// @brief Evaluates the final state invariant masses
    ///
    /// @tparam generator_t The random number generator type
    /// @param [in, out] generator The random number generator
    /// @param [in] kinematicParameters Parametrisation of kinematic properties
    ///
    /// @return Vector containing the invariant masses
    template <typename generator_t>
    std::vector<float> sampleInvariantMasses(generator_t& generator, const std::any& kinematicParmaeters) const;
  
    /// @brief Converter from sampled numbers to a vector of particles
    ///
    /// @param [in] pdgId The PDG IDs
    /// @param [in] momentum The momenta
    /// @param [in] invariantMass The invariant masses
    ///
    /// @return Vector containing the final state particles
	std::vector<Particle> convertParametersToParticles(const std::vector<int>& pdgId, const std::vector<float>& momentum, const std::vector<float>& invariantMass) const;
	
    ///Function gets random number rnd in the range [0,1) as argument 
    ///and returns function value according to a histogram distribution
    unsigned int sampleDiscreteValues(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const; // TODO: this should be only used for multiplicity
	
	
	double sampleContinuousValues(double rnd, const detail::Parameters::CumulativeDistribution& distribution) const;  
};

//~ vector<Event>
//~ runSimulation(const vector<Event>& events, const SimulationParameters& simParsLower, const SimulationParameters& simParsUpper, unsigned int mult, bool soft, unsigned int nSamples)
//~ {
	//~ vector<Event> result;
	//~ vector<pair<vector<double>, bool>> nonServedInvMass;
	//~ uniform_real_distribution<> d(0., 1.);

	//~ do
	//~ {
		//~ unsigned int samplesLeft = nonServedInvMass.empty() ? nSamples : nonServedInvMass.size();
//~ cout << "Current events: " << result.size() << endl;	
		//~ vector<vector<double>> momSamples, momSamplesTarget;
		//~ vector<vector<int>> pdgSamples;
		//~ vector<pair<vector<double>, bool>> invMassSamples;
		//~ vector<bool> dices;
		
		//~ for(unsigned int i = 0; i < samplesLeft; i++)
		//~ {
//~ if(i % 100 == 0)
//~ cout << i << " iteration" << endl;
			//~ double rnd = d(gen);
			//~ const bool dice = nonServedInvMass.empty() ? rnd < weightLower : nonServedInvMass[i].second;
			//~ dices.push_back(dice);
			//~ const SimulationParameters& simPars = dice ? simParsLower : simParsUpper;
			//~ auto sampleGauss = sampleInGauss(simPars.covAndMapper.first);
			//~ auto sampleMapped = mapBack(simPars.covAndMapper.second, sampleGauss, soft);
			//~ momSamples.push_back(sampleMapped);

			//~ vector<int> pdgSample;
			//~ for(unsigned int j = 0; j < sampleMapped.size(); j++)
				//~ if(j == 0)
				//~ {
					//~ pdgSample.push_back(samplePDG(simPars.pdgProb.at(events[0].pdg)));
				//~ }
				//~ else
				//~ {
					//~ pdgSample.push_back(samplePDG(simPars.pdgProb.at(pdgSample[j - 1])));
				//~ }
			//~ pdgSamples.push_back(pdgSample);
		//~ }
		//~ if(nonServedInvMass.empty())
		//~ {
			//~ for(unsigned int i = 0; i < samplesLeft; i++)
			//~ {
				//~ const SimulationParameters& simPars = dices[i] ? simParsLower : simParsUpper;
				//~ auto massGauss = sampleInGauss(simPars.invMassMapper.first);
				//~ auto massMapped = mapBackInvMass(simPars.invMassMapper.second, massGauss, soft);
				//~ invMassSamples.push_back(make_pair(massMapped, dices[i]));
			//~ }
		//~ }
		//~ else
		//~ {
			//~ invMassSamples = nonServedInvMass;
			//~ nonServedInvMass.clear();
		//~ }
		//~ momSamplesTarget = momSamples;
		//~ invertPostProcessMomDataToTaret(momSamplesTarget);
		//~ invertPostProcessMomData(momSamples, dices);
		//~ vector<Event> evtResult = convertSimulationToEvent(events, momSamples, momSamplesTarget, pdgSamples, invMassSamples, mult, soft, nonServedInvMass);
//~ cout << "Generated events: " << evtResult.size() << endl;
		//~ result.insert(result.end(), evtResult.begin(), evtResult.end());
	//~ } while(!nonServedInvMass.empty());
	//~ return result;
//~ }

}  // namespace ActsFatras
