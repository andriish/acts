// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include <cstdint>

class TH1;

class TFCS1DFunctionInt32Histogram
{
  public:
    TFCS1DFunctionInt32Histogram(const TH1* hist=nullptr);

    virtual void Initialize(const TH1* hist);
    
    typedef uint32_t HistoContent_t;
    static const HistoContent_t s_MaxValue;

    ///Function gets random number rnd in the range [0,1) as argument 
    ///and returns function value according to a histogram distribution
    virtual double rnd_to_fct(double rnd) const;

    const std::vector<float>& get_HistoBordersx() const {return m_HistoBorders;};
    std::vector<float>& get_HistoBordersx() {return m_HistoBorders;};
    const std::vector<HistoContent_t>& get_HistoContents() const {return m_HistoContents;};
    std::vector<HistoContent_t>& get_HistoContents() {return m_HistoContents;};
    
  private:
    
    std::vector<float> m_HistoBorders;
    std::vector<HistoContent_t> m_HistoContents;
};