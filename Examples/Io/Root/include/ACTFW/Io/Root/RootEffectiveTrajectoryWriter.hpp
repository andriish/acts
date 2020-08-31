// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/EventData/EffectiveMultiTrajectory.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

#include <mutex>
#include <vector>

class TFile;
class TTree;

namespace FW {

/// @class RootTrajectoryWriter
///
/// Write out a trajectory (i.e. a vector of
/// trackState at the moment) into a TTree
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// Each entry in the TTree corresponds to one trajectory for optimum
/// writing speed. The event number is part of the written data.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing
/// file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootEffectiveTrajectoryWriter final : public WriterT<std::vector<EffectiveMultiTrajectory>> {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::string inputParticles;     ///< input truth particles collection.
    std::string inputTrajectories;  ///< input (fitted) trajectories collection
    std::string outputDir;          ///< output directory
    std::string outputFilename = "tracks.root";  ///< output filename
    std::string outputTreename = "tracks";       ///< name of the output tree
    std::string fileMode = "RECREATE";           ///< file access mode
    TFile* rootFile = nullptr;                   ///< common root file
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  RootEffectiveTrajectoryWriter(const Config& cfg, Acts::Logging::Level lvl);

  /// Virtual destructor
  ~RootEffectiveTrajectoryWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trajectories are what to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const std::vector<EffectiveMultiTrajectory>& trajectories) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  int m_eventNr{0};              ///< the event number
  int m_trajNr{0};               ///< the trajectory number

  unsigned long m_t_barcode{0};  ///< Truth particle barcode
  int m_t_charge{0};             ///< Truth particle charge
  float m_t_time{0};             ///< Truth particle time
  float m_t_vx{-99.};            ///< Truth particle vertex x
  float m_t_vy{-99.};            ///< Truth particle vertex y
  float m_t_vz{-99.};            ///< Truth particle vertex z
  float m_t_px{-99.};            ///< Truth particle initial momentum px
  float m_t_py{-99.};            ///< Truth particle initial momentum py
  float m_t_pz{-99.};            ///< Truth particle initial momentum pz
  float m_t_theta{-99.};         ///< Truth particle initial momentum theta
  float m_t_phi{-99.};           ///< Truth particle initial momentum phi
  float m_t_pT{-99.};            ///< Truth particle initial momentum pT
  float m_t_eta{-99.};           ///< Truth particle initial momentum eta

  std::vector<float> m_t_x;  ///< Global truth hit position x
  std::vector<float> m_t_y;  ///< Global truth hit position y
  std::vector<float> m_t_z;  ///< Global truth hit position z
  std::vector<float> m_t_r;  ///< Global truth hit position r
  std::vector<float>
      m_t_dx;  ///< Truth particle direction x at global hit position
  std::vector<float>
      m_t_dy;  ///< Truth particle direction y at global hit position
  std::vector<float>
      m_t_dz;  ///< Truth particle direction z at global hit position

  std::vector<float> m_t_eLOC0;   ///< truth parameter eLOC_0
  std::vector<float> m_t_eLOC1;   ///< truth parameter eLOC_1
  std::vector<float> m_t_ePHI;    ///< truth parameter ePHI
  std::vector<float> m_t_eTHETA;  ///< truth parameter eTHETA
  std::vector<float> m_t_eQOP;    ///< truth parameter eQOP
  std::vector<float> m_t_eT;      ///< truth parameter eT

  int m_nStates{0};                 ///< number of all states
  int m_nMeasurements{0};           ///< number of states with measurements
  std::vector<int> m_volumeID;      ///< volume identifier
  std::vector<int> m_layerID;       ///< layer identifier
  std::vector<int> m_moduleID;      ///< surface identifier
  std::vector<float> m_lx_hit;      ///< uncalibrated measurement local x
  std::vector<float> m_ly_hit;      ///< uncalibrated measurement local y
  std::vector<float> m_x_hit;       ///< uncalibrated measurement global x
  std::vector<float> m_y_hit;       ///< uncalibrated measurement global y
  std::vector<float> m_z_hit;       ///< uncalibrated measurement global z
  std::vector<float> m_res_x_hit;   ///< hit residual x
  std::vector<float> m_res_y_hit;   ///< hit residual y
  std::vector<float> m_res_z_hit;   ///< hit residual y
  std::vector<float> m_err_x_hit;   ///< hit err x
  std::vector<float> m_err_y_hit;   ///< hit err y
  std::vector<float> m_err_z_hit;
  std::vector<float> m_pull_x_hit;  ///< hit pull x
  std::vector<float> m_pull_y_hit;  ///< hit pull y
  std::vector<float> m_pull_z_hit;  ///< hit pull y
  std::vector<int> m_dim_hit;       ///< dimension of measurement

  bool m_hasFittedParams;        ///< if the track has fitted parameter
  float m_eLOC0_fit{-99.};       ///< fitted parameter eLOC_0
  float m_eLOC1_fit{-99.};       ///< fitted parameter eLOC_1
  float m_ePHI_fit{-99.};        ///< fitted parameter ePHI
  float m_eTHETA_fit{-99.};      ///< fitted parameter eTHETA
  float m_eQOP_fit{-99.};        ///< fitted parameter eQOP
  float m_eT_fit{-99.};          ///< fitted parameter eT
  float m_err_eLOC0_fit{-99.};   ///< fitted parameter eLOC_-99.err
  float m_err_eLOC1_fit{-99.};   ///< fitted parameter eLOC_1 err
  float m_err_ePHI_fit{-99.};    ///< fitted parameter ePHI err
  float m_err_eTHETA_fit{-99.};  ///< fitted parameter eTHETA err
  float m_err_eQOP_fit{-99.};    ///< fitted parameter eQOP err
  float m_err_eT_fit{-99.};      ///< fitted parameter eT err

  int m_nPredicted{0};      ///< number of states with predicted parameter
  std::vector<bool> m_prt;  ///< predicted status
  std::vector<float> m_eLOC0_prt;       ///< predicted parameter eLOC0
  std::vector<float> m_eLOC1_prt;       ///< predicted parameter eLOC1
  std::vector<float> m_ePHI_prt;        ///< predicted parameter ePHI
  std::vector<float> m_eTHETA_prt;      ///< predicted parameter eTHETA
  std::vector<float> m_ePos0_prt;       ///< predicted parameter eLOC0
  std::vector<float> m_ePos1_prt;       ///< predicted parameter eLOC1
  std::vector<float> m_ePos2_prt;       ///< predicted parameter eLOC1
  std::vector<float> m_eDir0_prt;        ///< predicted parameter ePHI
  std::vector<float> m_eDir1_prt;      ///< predicted parameter eTHETA
  std::vector<float> m_eDir2_prt;      ///< predicted parameter eTHETA
  std::vector<float> m_eQOP_prt;        ///< predicted parameter eQOP
  std::vector<float> m_eT_prt;          ///< predicted parameter eT
  std::vector<float> m_res_eLOC0_prt;   ///< predicted parameter eLOC0 residual
  std::vector<float> m_res_eLOC1_prt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_res_ePHI_prt;    ///< predicted parameter ePHI residual
  std::vector<float> m_res_eTHETA_prt;  ///< predicted parameter eTHETA residual
  std::vector<float> m_res_ePos0_prt;   ///< predicted parameter eLOC0 residual
  std::vector<float> m_res_ePos1_prt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_res_ePos2_prt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_res_eDir0_prt;    ///< predicted parameter ePHI residual
  std::vector<float> m_res_eDir1_prt;  ///< predicted parameter eTHETA residual
  std::vector<float> m_res_eDir2_prt;  ///< predicted parameter eTHETA residual
  std::vector<float> m_res_eQOP_prt;    ///< predicted parameter eQOP residual
  std::vector<float> m_res_eT_prt;      ///< predicted parameter eT residual
  std::vector<float> m_err_eLOC0_prt;   ///< predicted parameter eLOC0 error
  std::vector<float> m_err_eLOC1_prt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_ePHI_prt;    ///< predicted parameter ePHI error
  std::vector<float> m_err_eTHETA_prt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_ePos0_prt;   ///< predicted parameter eLOC0 error
  std::vector<float> m_err_ePos1_prt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_ePos2_prt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_eDir0_prt;    ///< predicted parameter ePHI error
  std::vector<float> m_err_eDir1_prt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_eDir2_prt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_eQOP_prt;    ///< predicted parameter eQOP error
  std::vector<float> m_err_eT_prt;      ///< predicted parameter eT error
  std::vector<float> m_pull_eLOC0_prt;  ///< predicted parameter eLOC0 pull
  std::vector<float> m_pull_eLOC1_prt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_pull_ePHI_prt;   ///< predicted parameter ePHI pull
  std::vector<float> m_pull_eTHETA_prt;  ///< predicted parameter eTHETA pull
  std::vector<float> m_pull_ePos0_prt;  ///< predicted parameter eLOC0 pull
  std::vector<float> m_pull_ePos1_prt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_pull_ePos2_prt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_pull_eDir0_prt;   ///< predicted parameter ePHI pull
  std::vector<float> m_pull_eDir1_prt;  ///< predicted parameter eTHETA pull
  std::vector<float> m_pull_eDir2_prt;  ///< predicted parameter eTHETA pull
  std::vector<float> m_pull_eQOP_prt;    ///< predicted parameter eQOP pull
  std::vector<float> m_pull_eT_prt;      ///< predicted parameter eT pull
  std::vector<float> m_x_prt;            ///< predicted global x
  std::vector<float> m_y_prt;            ///< predicted global y
  std::vector<float> m_z_prt;            ///< predicted global z
  std::vector<float> m_px_prt;           ///< predicted momentum px
  std::vector<float> m_py_prt;           ///< predicted momentum py
  std::vector<float> m_pz_prt;           ///< predicted momentum pz
  std::vector<float> m_eta_prt;          ///< predicted momentum eta
  std::vector<float> m_pT_prt;           ///< predicted momentum pT

  int m_nFiltered{0};              ///< number of states with filtered parameter
  std::vector<bool> m_flt;  ///< predicted status
  std::vector<float> m_eLOC0_flt;       ///< predicted parameter eLOC0
  std::vector<float> m_eLOC1_flt;       ///< predicted parameter eLOC1
  std::vector<float> m_ePHI_flt;        ///< predicted parameter ePHI
  std::vector<float> m_eTHETA_flt;      ///< predicted parameter eTHETA
  std::vector<float> m_ePos0_flt;       ///< predicted parameter eLOC0
  std::vector<float> m_ePos1_flt;       ///< predicted parameter eLOC1
  std::vector<float> m_ePos2_flt;       ///< predicted parameter eLOC1
  std::vector<float> m_eDir0_flt;        ///< predicted parameter ePHI
  std::vector<float> m_eDir1_flt;      ///< predicted parameter eTHETA
  std::vector<float> m_eDir2_flt;      ///< predicted parameter eTHETA
  std::vector<float> m_eQOP_flt;        ///< predicted parameter eQOP
  std::vector<float> m_eT_flt;          ///< predicted parameter eT
  std::vector<float> m_res_eLOC0_flt;   ///< predicted parameter eLOC0 residual
  std::vector<float> m_res_eLOC1_flt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_res_ePHI_flt;    ///< predicted parameter ePHI residual
  std::vector<float> m_res_eTHETA_flt;  ///< predicted parameter eTHETA residual
  std::vector<float> m_res_ePos0_flt;   ///< predicted parameter eLOC0 residual
  std::vector<float> m_res_ePos1_flt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_res_ePos2_flt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_res_eDir0_flt;    ///< predicted parameter ePHI residual
  std::vector<float> m_res_eDir1_flt;  ///< predicted parameter eTHETA residual
  std::vector<float> m_res_eDir2_flt;  ///< predicted parameter eTHETA residual
  std::vector<float> m_res_eQOP_flt;    ///< predicted parameter eQOP residual
  std::vector<float> m_res_eT_flt;      ///< predicted parameter eT residual
  std::vector<float> m_err_eLOC0_flt;   ///< predicted parameter eLOC0 error
  std::vector<float> m_err_eLOC1_flt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_ePHI_flt;    ///< predicted parameter ePHI error
  std::vector<float> m_err_eTHETA_flt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_ePos0_flt;   ///< predicted parameter eLOC0 error
  std::vector<float> m_err_ePos1_flt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_ePos2_flt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_eDir0_flt;    ///< predicted parameter ePHI error
  std::vector<float> m_err_eDir1_flt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_eDir2_flt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_eQOP_flt;    ///< predicted parameter eQOP error
  std::vector<float> m_err_eT_flt;      ///< predicted parameter eT error
  std::vector<float> m_pull_eLOC0_flt;  ///< predicted parameter eLOC0 pull
  std::vector<float> m_pull_eLOC1_flt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_pull_ePHI_flt;   ///< predicted parameter ePHI pull
  std::vector<float> m_pull_eTHETA_flt;  ///< predicted parameter eTHETA pull
  std::vector<float> m_pull_ePos0_flt;  ///< predicted parameter eLOC0 pull
  std::vector<float> m_pull_ePos1_flt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_pull_ePos2_flt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_pull_eDir0_flt;   ///< predicted parameter ePHI pull
  std::vector<float> m_pull_eDir1_flt;  ///< predicted parameter eTHETA pull
  std::vector<float> m_pull_eDir2_flt;  ///< predicted parameter eTHETA pull
  std::vector<float> m_pull_eQOP_flt;    ///< predicted parameter eQOP pull
  std::vector<float> m_pull_eT_flt;      ///< predicted parameter eT pull
  std::vector<float> m_x_flt;            ///< predicted global x
  std::vector<float> m_y_flt;            ///< predicted global y
  std::vector<float> m_z_flt;            ///< predicted global z
  std::vector<float> m_px_flt;           ///< predicted momentum px
  std::vector<float> m_py_flt;           ///< predicted momentum py
  std::vector<float> m_pz_flt;           ///< predicted momentum pz
  std::vector<float> m_eta_flt;          ///< predicted momentum eta
  std::vector<float> m_pT_flt;           ///< predicted momentum pT
  std::vector<float> m_chi2;             ///< chisq from filtering

  int m_nSmoothed{0};              ///< number of states with smoothed parameter
  std::vector<bool> m_smt;  ///< predicted status
  std::vector<float> m_eLOC0_smt;       ///< predicted parameter eLOC0
  std::vector<float> m_eLOC1_smt;       ///< predicted parameter eLOC1
  std::vector<float> m_ePHI_smt;        ///< predicted parameter ePHI
  std::vector<float> m_eTHETA_smt;      ///< predicted parameter eTHETA
  std::vector<float> m_ePos0_smt;       ///< predicted parameter eLOC0
  std::vector<float> m_ePos1_smt;       ///< predicted parameter eLOC1
  std::vector<float> m_ePos2_smt;       ///< predicted parameter eLOC1
  std::vector<float> m_eDir0_smt;        ///< predicted parameter ePHI
  std::vector<float> m_eDir1_smt;      ///< predicted parameter eTHETA
  std::vector<float> m_eDir2_smt;      ///< predicted parameter eTHETA
  std::vector<float> m_eQOP_smt;        ///< predicted parameter eQOP
  std::vector<float> m_eT_smt;          ///< predicted parameter eT
  std::vector<float> m_res_eLOC0_smt;   ///< predicted parameter eLOC0 residual
  std::vector<float> m_res_eLOC1_smt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_res_ePHI_smt;    ///< predicted parameter ePHI residual
  std::vector<float> m_res_eTHETA_smt;  ///< predicted parameter eTHETA residual
  std::vector<float> m_res_ePos0_smt;   ///< predicted parameter eLOC0 residual
  std::vector<float> m_res_ePos1_smt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_res_ePos2_smt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_res_eDir0_smt;    ///< predicted parameter ePHI residual
  std::vector<float> m_res_eDir1_smt;  ///< predicted parameter eTHETA residual
  std::vector<float> m_res_eDir2_smt;  ///< predicted parameter eTHETA residual
  std::vector<float> m_res_eQOP_smt;    ///< predicted parameter eQOP residual
  std::vector<float> m_res_eT_smt;      ///< predicted parameter eT residual
  std::vector<float> m_err_eLOC0_smt;   ///< predicted parameter eLOC0 error
  std::vector<float> m_err_eLOC1_smt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_ePHI_smt;    ///< predicted parameter ePHI error
  std::vector<float> m_err_eTHETA_smt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_ePos0_smt;   ///< predicted parameter eLOC0 error
  std::vector<float> m_err_ePos1_smt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_ePos2_smt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_eDir0_smt;    ///< predicted parameter ePHI error
  std::vector<float> m_err_eDir1_smt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_eDir2_smt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_eQOP_smt;    ///< predicted parameter eQOP error
  std::vector<float> m_err_eT_smt;      ///< predicted parameter eT error
  std::vector<float> m_pull_eLOC0_smt;  ///< predicted parameter eLOC0 pull
  std::vector<float> m_pull_eLOC1_smt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_pull_ePHI_smt;   ///< predicted parameter ePHI pull
  std::vector<float> m_pull_eTHETA_smt;  ///< predicted parameter eTHETA pull
  std::vector<float> m_pull_ePos0_smt;  ///< predicted parameter eLOC0 pull
  std::vector<float> m_pull_ePos1_smt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_pull_ePos2_smt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_pull_eDir0_smt;   ///< predicted parameter ePHI pull
  std::vector<float> m_pull_eDir1_smt;  ///< predicted parameter eTHETA pull
  std::vector<float> m_pull_eDir2_smt;  ///< predicted parameter eTHETA pull
  std::vector<float> m_pull_eQOP_smt;    ///< predicted parameter eQOP pull
  std::vector<float> m_pull_eT_smt;      ///< predicted parameter eT pull
  std::vector<float> m_x_smt;            ///< predicted global x
  std::vector<float> m_y_smt;            ///< predicted global y
  std::vector<float> m_z_smt;            ///< predicted global z
  std::vector<float> m_px_smt;           ///< predicted momentum px
  std::vector<float> m_py_smt;           ///< predicted momentum py
  std::vector<float> m_pz_smt;           ///< predicted momentum pz
  std::vector<float> m_eta_smt;          ///< predicted momentum eta
  std::vector<float> m_pT_smt;           ///< predicted momentum pT
};

}  // namespace FW
