// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <stdexcept>
#include <string>

namespace FW {

/// Source link class for simulation in the acts-framework.
///
/// The source link stores the measuremts, surface, and the associated simulated
/// truth hit.
///
/// @todo Allow multiple truth hits e.g. for merged hits.
class EffectiveSourceLink {
 public:
  EffectiveSourceLink(const Acts::GeometryObject& referenceObject, const ActsFatras::Hit& truthHit,
                size_t dim, Acts::BoundVector values, Acts::BoundMatrix cov)
      : m_values(values), m_cov(cov), m_dim(dim),
        m_geometryId(truthHit.geometryId()),
        m_referenceObject(&referenceObject),
        m_truthHit(&truthHit) {}
  
  /// Must be default_constructible to satisfy SourceLinkConcept.
  EffectiveSourceLink() = default;
  EffectiveSourceLink(EffectiveSourceLink&&) = default;
  EffectiveSourceLink(const EffectiveSourceLink&) = default;
  EffectiveSourceLink& operator=(EffectiveSourceLink&&) = default;
  EffectiveSourceLink& operator=(const EffectiveSourceLink&) = default;
  ~EffectiveSourceLink() { if(m_meas != nullptr) delete(m_meas); }

  using MeasurementType = std::variant<Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices, Acts::BoundParametersIndices::eLOC_0>,
						  Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,Acts::BoundParametersIndices::eLOC_0, Acts::BoundParametersIndices::eLOC_1>, 
						  Acts::Measurement<EffectiveSourceLink, Acts::FreeParametersIndices, Acts::FreeParametersIndices::eFreePos0, Acts::FreeParametersIndices::eFreePos1, Acts::FreeParametersIndices::eFreePos2>>;
						  
  constexpr Acts::GeometryID geometryId() const { return m_geometryId; }
  constexpr const Acts::GeometryObject& referenceObject() const { return *m_referenceObject; }
  constexpr const ActsFatras::Hit& truthHit() const { return *m_truthHit; }
  MeasurementType operator*() const 
  { 
	  	if (m_dim == 0) {
		  throw std::runtime_error("Cannot create dim 0 measurement");
		} else if (m_dim == 1) {
		  return Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,
								   Acts::ParDef::eLOC_0>(dynamic_cast<const Acts::Surface*>(m_referenceObject)->getSharedPtr(), *this, m_cov.topLeftCorner<1, 1>(), m_values[0]);
		} else if (m_dim == 2) {
		  return Acts::Measurement<EffectiveSourceLink, Acts::BoundParametersIndices,
								   Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_1>(
			  dynamic_cast<const Acts::Surface*>(m_referenceObject)->getSharedPtr(), *this, m_cov.topLeftCorner<2, 2>(),
			  m_values[0], m_values[1]);
	    } else if (m_dim == 3) {
			return Acts::Measurement<EffectiveSourceLink, Acts::FreeParametersIndices,
								   Acts::FreeParametersIndices::eFreePos0, Acts::FreeParametersIndices::eFreePos1, Acts::FreeParametersIndices::eFreePos2>(
			  dynamic_cast<const Acts::Volume*>(m_referenceObject)->getSharedPtr(), *this, m_cov.topLeftCorner<3, 3>(),
			  m_values[0], m_values[1], m_values[2]);
		} else {
		  throw std::runtime_error("Dim " + std::to_string(m_dim) +
								   " currently not supported.");
		}
	}
	  
 private:
   Acts::BoundVector m_values;
  Acts::BoundMatrix m_cov;
  size_t m_dim = 0u;
	 Acts::GeometryID m_geometryId;
	 const Acts::GeometryObject* m_referenceObject;
	 MeasurementType* m_meas;
	 const ActsFatras::Hit* m_truthHit;
  

  friend constexpr bool operator==(const EffectiveSourceLink& lhs,
                                   const EffectiveSourceLink& rhs) {
    return lhs.m_truthHit == rhs.m_truthHit;
  }
};

/// Store source links ordered by geometry identifier.
using EffectiveSourceLinkContainer = GeometryIdMultiset<EffectiveSourceLink>;

}  // end of namespace FW