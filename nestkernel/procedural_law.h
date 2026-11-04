/*
 *  procedural_law.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*
 * EXPERIMENTAL, NOT STOCK NEST (regenerated connections, Stage A).
 *
 * Resolves a projection's weight or delay, as the connection builder holds it, into a ProceduralLaw
 * (procedural_random.h), once at building time. Delivery then draws the law directly instead of calling NEST's
 * parameter objects, which costs 2.3 times as much (stage_a/README.md in the scoping repository).
 *
 * Accepted: a scalar, or a constant, normal or lognormal parameter, alone or inside one math.max, math.min or
 * math.redraw. These are the forms the macaque and human multi-area models use; anything else is refused.
 */

#ifndef PROCEDURAL_LAW_H
#define PROCEDURAL_LAW_H

// C++ includes:
#include <memory>

// Includes from nestkernel:
#include "conn_parameter.h"
#include "exceptions.h"
#include "parameter.h"
#include "procedural_random.h"

namespace nest
{

inline ProceduralLaw
resolve_procedural_law( const ConnParameter& cp )
{
  ProceduralLaw law;
  if ( cp.is_scalar() )
  {
    law.mean = cp.value_double( 0, nullptr, 0, nullptr );
    return law;
  }
  const auto* wrapper = dynamic_cast< const ParameterConnParameterWrapper* >( &cp );
  if ( wrapper == nullptr )
  {
    throw BadParameter( "Regenerated connections take a scalar or a parameter object as weight and delay." );
  }

  ParameterPTR p = wrapper->get_parameter();
  if ( const auto mx = std::dynamic_pointer_cast< MaxParameter >( p ) )
  {
    law.clip = ProceduralLaw::Clip::max;
    law.low = mx->get_bound();
    p = mx->get_operand();
  }
  else if ( const auto mn = std::dynamic_pointer_cast< MinParameter >( p ) )
  {
    law.clip = ProceduralLaw::Clip::min;
    law.high = mn->get_bound();
    p = mn->get_operand();
  }
  else if ( const auto rd = std::dynamic_pointer_cast< RedrawParameter >( p ) )
  {
    law.clip = ProceduralLaw::Clip::redraw;
    law.low = rd->get_min();
    law.high = rd->get_max();
    law.max_redraws = rd->get_max_redraws();
    p = rd->get_operand();
  }

  if ( const auto c = std::dynamic_pointer_cast< ConstantParameter >( p ) )
  {
    law.base = ProceduralLaw::Base::constant;
    law.mean = c->get_value();
  }
  else if ( const auto n = std::dynamic_pointer_cast< NormalParameter >( p ) )
  {
    law.base = ProceduralLaw::Base::normal;
    law.mean = n->get_mean();
    law.sigma = n->get_std();
  }
  else if ( const auto ln = std::dynamic_pointer_cast< LognormalParameter >( p ) )
  {
    law.base = ProceduralLaw::Base::lognormal;
    law.mean = ln->get_mean();
    law.sigma = ln->get_std();
  }
  else
  {
    throw BadParameter(
      "Regenerated connections take a constant, normal or lognormal weight or delay, alone or inside one math.max, "
      "math.min or math.redraw." );
  }
  return law;
}

} // namespace nest

#endif /* #ifndef PROCEDURAL_LAW_H */
