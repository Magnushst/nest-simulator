/*
 *  dict_util.h
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

#ifndef DICT_UTIL_H
#define DICT_UTIL_H

// C++ includes:
#include <cstdint>
#include <random>

// Includes from nestkernel:
#include "kernel_manager.h"
#include "random_generators.h" // EXPERIMENTAL, NOT STOCK NEST (the invariant branch: RandomGenerator)
#include "vp_manager_impl.h"

#include "dictionary.h"

namespace nest
{
/**
 * Obtain value from parameter dictionary including evaluation of random or spatial parameters.
 *
 * This function should be used instead of `update_value()` everywhere where the user may pass
 * random or spatial parameters.
 */
template < typename T >
bool
update_value_param( Dictionary const& d, const std::string& key, T& value, nest::Node* node )
{
  assert( node != nullptr );
  const auto it = d.find( key );
  if ( it != d.end() and std::holds_alternative< std::shared_ptr< nest::Parameter > >( it->second.item ) )
  {
    const auto param = d.get< ParameterPTR >( key );
    // EXPERIMENTAL, NOT STOCK NEST (the invariant branch, NEST_INVARIANT.md): drawn from a generator keyed by the base
    // seed, the node and the parameter's name (FNV-1a), not by the node's virtual process, so a node gets the same
    // value at any number of processes
    std::uint32_t name_hash = 2166136261u;
    for ( const char c : key )
    {
      name_hash = ( name_hash ^ static_cast< unsigned char >( c ) ) * 16777619u;
    }
    const std::uint64_t id = node->get_node_id();
    RandomGenerator< std::mt19937_64 > keyed( { static_cast< std::uint32_t >( kernel().random_manager.get_base_seed() ),
      static_cast< std::uint32_t >( id ),
      static_cast< std::uint32_t >( id >> 32 ),
      name_hash } );
    value = param->value( &keyed, node );
    return true;
  }
  else
  {
    return d.update_value( key, value );
  }
}

}  // namespace nest

#endif  // DICT_UTIL_H
