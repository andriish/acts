// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootNuclearInteractionParametersReader.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iostream>

#include <TFile.h>
#include <TH1F.h>
#include <TVectorF.h>
#include <ActsFatras/Physics/NuclearInteraction/Parameters.hpp>

namespace {
	
/// @brief This method builds components for transforming a probability distribution into components that represent the cumulative probability distribution
///
/// @param [in] hist The probability distribution
///
/// @return Tuple containing the bin borders, the non-normalised cumulative distribution and the sum over all entries
std::tuple<std::vector<float>, std::vector<double>, double>
buildNotNormalisedMap(TH1F const* hist) {
  // Retrieve the number of bins & borders
  const int nBins = hist->GetNbinsX();
  std::vector<float> histoBorders(nBins + 1);

  // Fill the cumulative histogram
  double integral = 0.;
  std::vector<double> temp_HistoContents(nBins);
  int iBin;
  for (iBin = 0; iBin < nBins; iBin++) {
    float binval = hist->GetBinContent(iBin + 1);
    // Avoid negative bin values
    if (binval < 0) {
      binval = 0.;
    }
    // Store the value
    integral += binval;
    temp_HistoContents[iBin] = integral;
  }

  // Ensure that content is available
  if (integral == 0.) {
    histoBorders.clear();
    temp_HistoContents.clear();
    return std::make_tuple(histoBorders, temp_HistoContents, integral);
  }

  // Set the bin borders
  for (iBin = 1; iBin <= nBins; iBin++)
    histoBorders[iBin - 1] = hist->GetXaxis()->GetBinLowEdge(iBin);
  histoBorders[nBins] = hist->GetXaxis()->GetXmax();

  return std::make_tuple(histoBorders, temp_HistoContents, integral);
}

/// @brief This method transforms a probability distribution into components
/// that represent the cumulative probability distribution
///
/// @note This method is used to avoid Root dependencies afterwards by
/// decomposing the histogram
/// @param [in] hist The probability distribution
///
/// @return Pair containing the bin borders and the bin content
std::pair<std::vector<float>, std::vector<uint32_t>> buildMap(TH1F const* hist) {
  // Build the components
  std::tuple<std::vector<float>, std::vector<double>, double> map = buildNotNormalisedMap(hist);
  const int nBins = hist->GetNbinsX();
  std::vector<double>& histoContents = std::get<1>(map);
  
  // Fast exit if the histogram is empty
  if(histoContents.empty())
	return std::make_pair(std::get<0>(map), std::vector<uint32_t>());

  // Set the bin content
  std::vector<uint32_t> normalisedHistoContents(nBins);
  const double invIntegral = 1. / std::get<2>(map);
  for (int iBin = 0; iBin < nBins; ++iBin) {
    normalisedHistoContents[iBin] = UINT32_MAX * (histoContents[iBin] * invIntegral);
  }
  
  return std::make_pair(std::get<0>(map), normalisedHistoContents);
}

/// @brief This method transforms a probability distribution into components
/// that represent the cumulative probability distribution
///
/// @note This method is used to avoid Root dependencies afterwards by
/// decomposing the histogram
/// @param [in] hist The probability distribution
/// @param [in] integral Scaling factor of the distribution
///
/// @return Pair containing the bin borders and the bin content
std::pair<std::vector<float>, std::vector<uint32_t>> buildMap(TH1F const* hist, double integral) {
  // Build the components
  std::tuple<std::vector<float>, std::vector<double>, double> map = buildNotNormalisedMap(hist);
  const int nBins = hist->GetNbinsX();
  std::vector<double>& histoContents = std::get<1>(map);
  
  // Fast exit if the histogram is empty
  if(histoContents.empty())
	return std::make_pair(std::get<0>(map), std::vector<uint32_t>());

  // Set the bin content
  std::vector<uint32_t> normalisedHistoContents(nBins);
  const double invIntegral = 1. / integral;
  for (int iBin = 0; iBin < nBins; ++iBin) {
    normalisedHistoContents[iBin] = UINT32_MAX * (histoContents[iBin] * invIntegral);
  }

  return std::make_pair(std::get<0>(map), normalisedHistoContents);
}

}
ActsExamples::RootNuclearInteractionParametersReader::RootNuclearInteractionParametersReader(
    const ActsExamples::RootNuclearInteractionParametersReader::Config& cfg)
    : ActsExamples::IReader(), m_cfg(cfg) {

  //~ // loop over the input files
  //~ for (auto inputFile : m_cfg.fileList) {
    //~ // add file to the input chain
    //~ m_inputChain->Add(inputFile.c_str());
    //~ ACTS_DEBUG("Adding File " << inputFile << " to tree '" << m_cfg.treeName
                              //~ << "'.");
  //~ }

  //~ m_events = m_inputChain->GetEntries();
  //~ ACTS_DEBUG("The full chain has " << m_events << " entries.");
}

ActsExamples::RootNuclearInteractionParametersReader::~RootNuclearInteractionParametersReader() {
}

std::string ActsExamples::RootNuclearInteractionParametersReader::name() const {
  return m_cfg.name;
}

std::pair<size_t, size_t>
ActsExamples::RootNuclearInteractionParametersReader::availableEvents() const {
  return {0u, 1u};
}

std::vector<std::pair<std::vector<float>, std::vector<uint32_t>>> ActsExamples::RootNuclearInteractionParametersReader::buildMaps(
    const std::vector<TH1F*>& histos) const {
  std::vector<std::pair<std::vector<float>, std::vector<uint32_t>>> maps;
  for (auto& h : histos) {
    maps.push_back(buildMap(h));
  }
  return maps;
}


//~ vector<SimulationParameters>
//~ loadSimulationParameters(unsigned int multMin, unsigned int multMax, string direc)
//~ {
	//~ vector<SimulationParameters> result;

	//~ TFile tf((direc + "simPars.root").c_str(), "read");
	//~ for(unsigned int i = multMin; i < multMax; i++)
	//~ {
		//~ unsigned int upToInvMass = i;
		//~ {
			//~ string nameSoft = "Soft_";
			//~ nameSoft += to_string(i);
					
			//~ unordered_map<int, unordered_map<int, double>> pdgProb;
			//~ unordered_map<int, double>* hist;
			//~ for(const auto& p1 : pdgToMass)
			//~ {
				//~ gDirectory->GetObject(("PDGprobability_" + to_string(p1.first) + "_" + nameSoft).c_str(), hist);
				//~ pdgProb[p1.first] = *hist;
			//~ }
		
			//~ TVectorD* eValues;
			//~ TMatrixD* eVectors;
			//~ TVectorD* rotParMeans;
			//~ gDirectory->GetObject(("Covariance_Momentum_Eigenvalues_" + nameSoft).c_str(), eValues);
			//~ gDirectory->GetObject(("Covariance_Momentum_Eigenvectors_" + nameSoft).c_str(), eVectors);
			//~ gDirectory->GetObject(("Covariance_Momentum_rotParMeans_" + nameSoft).c_str(), rotParMeans);
			//~ Covariance* cov = new Covariance(*eValues, *eVectors, *rotParMeans); // Memory leak?

			//~ unordered_map<string, TFCS1DFunctionInt32Histogram> momMapper;
			//~ string basenameSoft = "MomentumSoft_";
			//~ basenameSoft += to_string(i) + "_";
			//~ for(unsigned int j = 0; j < i + 1; j++)
			//~ {
				//~ string nameSoftHist = basenameSoft + to_string(j);
				//~ TFCS1DFunctionInt32Histogram* hist;
				//~ gDirectory->GetObject(("MomentumMaps_" + nameSoftHist + "_" + nameSoft).c_str(), hist);
				//~ momMapper[nameSoftHist] = *hist;
			//~ }
			//~ auto covAndMapper = make_pair(*cov, momMapper);
			
			//~ gDirectory->GetObject(("Covariance_InvariantMass_Eigenvalues_" + nameSoft).c_str(), eValues);
			//~ gDirectory->GetObject(("Covariance_InvariantMass_Eigenvectors_" + nameSoft).c_str(), eVectors);
			//~ gDirectory->GetObject(("Covariance_InvariantMass_rotParMeans_" + nameSoft).c_str(), rotParMeans);
			//~ cov = new Covariance(*eValues, *eVectors, *rotParMeans); // Memory leak?
			
			//~ vector<TFCS1DFunctionInt32Histogram> invMassMapper;
			//~ for(unsigned int j = 0; j < upToInvMass; j++)
			//~ {
				//~ string nameS = "InvariantMass_Maps_" + to_string(j) + "_Soft_" + to_string(i);
				//~ TFCS1DFunctionInt32Histogram* hist;
				//~ gDirectory->GetObject((nameS).c_str(), hist);
				//~ invMassMapper.push_back(*hist);
			//~ }
			//~ auto invariantMassMapper = make_pair(*cov, invMassMapper);
			
			//~ result.emplace_back(pdgProb, covAndMapper, invariantMassMapper, i, true);
		//~ }
		//~ {
			//~ string nameHard = to_string(i);
			
			//~ unordered_map<int, unordered_map<int, double>> pdgProb;
			//~ unordered_map<int, double>* hist;
			//~ for(const auto& p1 : pdgToMass)
			//~ {
				//~ gDirectory->GetObject(("PDGprobability_" + to_string(p1.first) + "_" + nameHard).c_str(), hist);
				//~ pdgProb[p1.first] = *hist;
			//~ }
		
			//~ TVectorD* eValues;
			//~ TMatrixD* eVectors;
			//~ TVectorD* rotParMeans;
			//~ gDirectory->GetObject(("Covariance_Momentum_Eigenvalues_" + nameHard).c_str(), eValues);
			//~ gDirectory->GetObject(("Covariance_Momentum_Eigenvectors_" + nameHard).c_str(), eVectors);
			//~ gDirectory->GetObject(("Covariance_Momentum_rotParMeans_" + nameHard).c_str(), rotParMeans);
			//~ Covariance* cov = new Covariance(*eValues, *eVectors, *rotParMeans);

			//~ unordered_map<string, TFCS1DFunctionInt32Histogram> momMapper;
			//~ string basenameHard = "Momentum_";
			//~ basenameHard += to_string(i) + "_";
			//~ for(unsigned int j = 0; j < i + 1; j++)
			//~ {
				//~ string nameHardHist = basenameHard + to_string(j);
				//~ TFCS1DFunctionInt32Histogram* hist;
				//~ gDirectory->GetObject(("MomentumMaps_" + nameHardHist + "_" + nameHard).c_str(), hist);
				//~ momMapper[nameHardHist] = *hist;
			//~ }
			//~ auto covAndMapper = make_pair(*cov, momMapper);
			
			//~ gDirectory->GetObject(("Covariance_InvariantMass_Eigenvalues_" + nameHard).c_str(), eValues);
			//~ gDirectory->GetObject(("Covariance_InvariantMass_Eigenvectors_" + nameHard).c_str(), eVectors);
			//~ gDirectory->GetObject(("Covariance_InvariantMass_rotParMeans_" + nameHard).c_str(), rotParMeans);
			//~ cov = new Covariance(*eValues, *eVectors, *rotParMeans); // Memory leak?
			
			//~ vector<TFCS1DFunctionInt32Histogram> invMassMapper;
			//~ for(unsigned int j = 0; j < upToInvMass; j++)
			//~ {
				//~ string nameH = "InvariantMass_Maps_" + to_string(j) + "_" + to_string(i);
				//~ TFCS1DFunctionInt32Histogram* hist;
				//~ gDirectory->GetObject((nameH).c_str(), hist);
				//~ invMassMapper.push_back(*hist);
			//~ }
			//~ auto invariantMassMapper = make_pair(*cov, invMassMapper);
			
			//~ result.emplace_back(pdgProb, covAndMapper, invariantMassMapper, i, true);
		//~ }
	//~ }
	//~ tf.Close();
	//~ return result;
//~ }


ActsExamples::ProcessCode ActsExamples::RootNuclearInteractionParametersReader::read(
    const ActsExamples::AlgorithmContext& context) {
	
	ACTS_DEBUG("Trying to read nulcear interaction parametrisations.");
	
  // Read and prepare the parametrisation
  if (context.eventNumber < 1) {
    // lock the mutex
    std::lock_guard<std::mutex> lock(m_read_mutex);
	
    // now read
	for(const std::string& file : m_cfg.fileList)
	{
		TFile tf(file.c_str(), "read");
		gDirectory->cd();
		auto list = gDirectory->GetListOfKeys();
		auto elem = list->First();
		while(elem)
		{
			ActsFatras::detail::Parameters parameters;
			
			char const* name = elem->GetName();
			parameters.momentum = std::stof(name);
			gDirectory->cd(name);

			TH1F* nuclearInteraction = (TH1F*) gDirectory->Get("NuclearInteraction");
			parameters.nuclearInteractionProbability = buildMap(nuclearInteraction, m_cfg.nSimulatedEvents);

			TVectorF* softInteraction = (TVectorF*) gDirectory->Get("SoftInteraction");
			parameters.softInteractionProbability = (*softInteraction)[0];

			std::vector<int> branchingPdgIds = *((std::vector<int>*) gDirectory->Get("BranchingPdgIds"));
			std::vector<int> targetPdgIds = *((std::vector<int>*) gDirectory->Get("TargetPdgIds"));
			std::vector<float> targetPdgProbability = *((std::vector<float>*) gDirectory->Get("TargetPdgProbability"));
			for(unsigned int i = 0; i < branchingPdgIds.size(); i++)
				parameters.pdgMap[branchingPdgIds[i]][targetPdgIds[i]] = targetPdgProbability[i];


			gDirectory->cd("soft");
			TH1F* softMultiplicity = (TH1F*) gDirectory->Get("Multiplicity");
			parameters.softMultiplicity = buildMap(softMultiplicity);

			auto softList = gDirectory->GetListOfKeys();
			auto softElement = softList->First();
			while(softElement)
			{
				if(softElement->IsFolder())
				{
					const char* distributionName = softElement->GetName();
					unsigned int mult = std::stoi(distributionName);
					gDirectory->cd(distributionName);
					std::vector<TH1F*> momentumDistributions;
					momentumDistributions.reserve(mult + 1);
					std::vector<TH1F*> invariantMassDistributions;
					invariantMassDistributions.reserve(mult);
					for(unsigned int i = 0; i < mult; i++)
					{
						momentumDistributions.push_back(std::move((TH1F*) gDirectory->Get(("MomentumDistribution_" + std::to_string(i)).c_str())));
						invariantMassDistributions.push_back(std::move((TH1F*) gDirectory->Get(("InvariantMassDistribution_" + std::to_string(i)).c_str())));
					}
					momentumDistributions.push_back(std::move((TH1F*) gDirectory->Get(("MomentumDistribution_" + std::to_string(mult)).c_str())));
					
					if(mult >= parameters.softMomentumDistributions.size())
						parameters.softMomentumDistributions.resize(mult + 1);
					if(mult >= parameters.softInvariantMassDistributions.size())
						parameters.softInvariantMassDistributions.resize(mult + 1);
					
					parameters.softMomentumDistributions[mult] = buildMaps(momentumDistributions);
					parameters.softInvariantMassDistributions[mult] = buildMaps(invariantMassDistributions);

					gDirectory->cd("..");
				}
				softElement = softList->After(softElement);
			}
			
			
			gDirectory->cd("../hard");
			TH1F* hardMultiplicity = (TH1F*) gDirectory->Get("Multiplicity");
			parameters.hardMultiplicity = buildMap(hardMultiplicity);

			auto hardList = gDirectory->GetListOfKeys();
			auto hardElement = hardList->First();
			while(hardElement)
			{
				if(hardElement->IsFolder())
				{
					const char* distributionName = hardElement->GetName();
					unsigned int mult = std::stoi(distributionName);
					gDirectory->cd(distributionName);
					std::vector<TH1F*> momentumDistributions;
					momentumDistributions.reserve(mult + 1);
					std::vector<TH1F*> invariantMassDistributions;
					invariantMassDistributions.reserve(mult);
					for(unsigned int i = 0; i < mult; i++)
					{
						momentumDistributions.push_back(std::move((TH1F*) gDirectory->Get(("MomentumDistribution_" + std::to_string(i)).c_str())));
						invariantMassDistributions.push_back(std::move((TH1F*) gDirectory->Get(("InvariantMassDistribution_" + std::to_string(i)).c_str())));
					}
					momentumDistributions.push_back(std::move((TH1F*) gDirectory->Get(("MomentumDistribution_" + std::to_string(mult)).c_str())));
					
					if(mult >= parameters.hardMomentumDistributions.size())
						parameters.hardMomentumDistributions.resize(mult + 1);
					if(mult >= parameters.hardInvariantMassDistributions.size())
						parameters.hardInvariantMassDistributions.resize(mult + 1);
					
					parameters.hardMomentumDistributions[mult] = buildMaps(momentumDistributions);
					parameters.hardInvariantMassDistributions[mult] = buildMaps(invariantMassDistributions);

					gDirectory->cd("..");
				}
				hardElement = softList->After(hardElement);
			}
			elem = list->After(elem); // TODO: this might be not needed
		}
	}
		
  }
  
  //~ // read in the material track
  //~ if (m_inputChain && context.eventNumber < m_events) {
    //~ // lock the mutex
    //~ std::lock_guard<std::mutex> lock(m_read_mutex);
    //~ // now read

    //~ // The collection to be written
    //~ std::vector<Acts::RecordedMaterialTrack> mtrackCollection;

    //~ for (size_t ib = 0; ib < m_cfg.batchSize; ++ib) {
      //~ // Read the correct entry: batch size * event_number + ib
      //~ m_inputChain->GetEntry(m_cfg.batchSize * context.eventNumber + ib);
      //~ ACTS_VERBOSE("Reading entry: " << m_cfg.batchSize * context.eventNumber +
                                            //~ ib);

      //~ Acts::RecordedMaterialTrack rmTrack;
      //~ // Fill the position and momentum
      //~ rmTrack.first.first = Acts::Vector3D(m_v_x, m_v_y, m_v_z);
      //~ rmTrack.first.second = Acts::Vector3D(m_v_px, m_v_py, m_v_pz);

      //~ // Fill the individual steps
      //~ size_t msteps = m_step_length->size();
      //~ ACTS_VERBOSE("Reading " << msteps << " material steps.");
      //~ rmTrack.second.materialInteractions.reserve(msteps);
      //~ rmTrack.second.materialInX0 = 0.;
      //~ rmTrack.second.materialInL0 = 0.;

      //~ for (size_t is = 0; is < msteps; ++is) {
        //~ double mX0 = (*m_step_X0)[is];
        //~ double mL0 = (*m_step_L0)[is];
        //~ double s = (*m_step_length)[is];

        //~ rmTrack.second.materialInX0 += s / mX0;
        //~ rmTrack.second.materialInL0 += s / mL0;

        //~ /// Fill the position & the material
        //~ Acts::MaterialInteraction mInteraction;
        //~ mInteraction.position =
            //~ Acts::Vector3D((*m_step_x)[is], (*m_step_y)[is], (*m_step_z)[is]);
        //~ mInteraction.direction = Acts::Vector3D(
            //~ (*m_step_dx)[is], (*m_step_dy)[is], (*m_step_dz)[is]);
        //~ mInteraction.materialSlab = Acts::MaterialSlab(
            //~ Acts::Material::fromMassDensity(mX0, mL0, (*m_step_A)[is],
                                            //~ (*m_step_Z)[is], (*m_step_rho)[is]),
            //~ s);
        //~ rmTrack.second.materialInteractions.push_back(std::move(mInteraction));
      //~ }
      //~ mtrackCollection.push_back(std::move(rmTrack));
    //~ }

    //~ // Write to the collection to the EventStore
    //~ context.eventStore.add(m_cfg.outputParametrisationStem, std::move(mtrackCollection));
  //~ }
  // Return success flag
  return ActsExamples::ProcessCode::SUCCESS;
}
