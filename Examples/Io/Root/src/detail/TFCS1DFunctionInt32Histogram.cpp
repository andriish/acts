// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TFCS1DFunctionInt32Histogram.hpp"
//~ #include <algorithm>
#include "TMath.h"
#include "TH1.h"

#include <iostream>

const TFCS1DFunctionInt32Histogram::HistoContent_t TFCS1DFunctionInt32Histogram::s_MaxValue=UINT32_MAX;

TFCS1DFunctionInt32Histogram::TFCS1DFunctionInt32Histogram(const TH1* hist=nullptr) 
{if(hist) Initialize(hist);}

void TFCS1DFunctionInt32Histogram::Initialize(const TH1* hist)
{
  Int_t nbinsx = hist->GetNbinsX();
  Int_t nbins  = nbinsx;
  
  float integral=0;
  m_HistoBorders.resize(nbinsx+1);
  m_HistoContents.resize(nbins);
  std::vector<double> temp_HistoContents(nbins);  
  int ibin=0;
  for (int ix=1; ix<=nbinsx; ix++){
    float binval=hist->GetBinContent(ix);
    if(binval<0) {
      //Can't work if a bin is negative, forcing bins to 0 in this case
      double fraction=binval/hist->Integral();
      if(TMath::Abs(fraction)>1e-5) {
        std::cout<<"WARNING: bin content is negative in histogram "<<hist->GetName()<<" : "<<hist->GetTitle()<<" binval="<<binval<<" "<<fraction*100<<"% of integral="<<hist->Integral()<<". Forcing bin to 0."<<std::endl;
      }  
      binval=0;
    }
    integral+=binval;
    temp_HistoContents[ibin]=integral;
    ++ibin;
  }
  if(integral<=0) {
    std::cout<<"ERROR: histogram "<<hist->GetName()<<" : "<<hist->GetTitle()<<" integral="<<integral<<" is <=0"<<std::endl;
    m_HistoBorders.resize(0);
    m_HistoContents.resize(0);
    return;
  }

  for (int ix=1; ix<=nbinsx; ix++) m_HistoBorders[ix-1]=hist->GetXaxis()->GetBinLowEdge(ix);
  m_HistoBorders[nbinsx]=hist->GetXaxis()->GetXmax();
  
  for(ibin=0;ibin<nbins;++ibin) {
    m_HistoContents[ibin]=s_MaxValue*(temp_HistoContents[ibin]/integral);
  }
}

double TFCS1DFunctionInt32Histogram::rnd_to_fct(double rnd) const
{
  if(m_HistoContents.size()==0) {
    return 0;
  }
  HistoContent_t int_rnd=s_MaxValue*rnd;
  auto it = std::upper_bound(m_HistoContents.begin(),m_HistoContents.end(),int_rnd);
  int ibin=std::distance(m_HistoContents.begin(),it);
  if(ibin>=(int)m_HistoContents.size()) ibin=m_HistoContents.size()-1;
  Int_t binx = ibin;
  
  HistoContent_t basecont=0;
  if(ibin>0) basecont=m_HistoContents[ibin-1];
  
  HistoContent_t dcont=m_HistoContents[ibin]-basecont;
  if(dcont>0) {
    return m_HistoBorders[binx] + ((m_HistoBorders[binx+1]-m_HistoBorders[binx]) * (int_rnd-basecont)) / dcont;
  } else {                             
    return m_HistoBorders[binx] + (m_HistoBorders[binx+1]-m_HistoBorders[binx]) / 2;
  }
}