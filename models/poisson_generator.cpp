/*
 *  poisson_generator.cpp
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

#include "poisson_generator.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "order_certified_drive.h"  // EXPERIMENTAL-NOT-STOCK (order-certified step)
#include "nest_impl.h"

// Includes from libnestutil:
#include "dict_util.h"

void
nest::register_poisson_generator( const std::string& name )
{
  register_node_model< poisson_generator >( name );
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameter
 * ---------------------------------------------------------------- */

nest::poisson_generator::Parameters_::Parameters_()
  : rate_( 0.0 )  // spks/s
{
}


/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::poisson_generator::Parameters_::get( Dictionary& d ) const
{
  d[ names::rate ] = rate_;
}

void
nest::poisson_generator::Parameters_::set( const Dictionary& d, Node* node )
{
  update_value_param( d, names::rate, rate_, node );
  if ( rate_ < 0 )
  {
    throw BadProperty( "The rate cannot be negative." );
  }
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::poisson_generator::poisson_generator()
  : StimulationDevice()
  , P_()
{
}

nest::poisson_generator::poisson_generator( const poisson_generator& n )
  : StimulationDevice( n )
  , P_( n.P_ )
{
}


/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::poisson_generator::init_state_()
{
  StimulationDevice::init_state();
}

void
nest::poisson_generator::init_buffers_()
{
  StimulationDevice::init_buffers();
}

void
nest::poisson_generator::pre_run_hook()
{
  StimulationDevice::pre_run_hook();

  // rate_ is in spks/s, dt in ms, so we have to convert from s to ms
  poisson_distribution::param_type param( Time::get_resolution().get_ms() * P_.rate_ * 1e-3 );
  V_.poisson_dist_.param( param );

  // EXPERIMENTAL-NOT-STOCK (order-certified step)
  order_certified_drive::register_hook();
  order_certified_drive::mark_dirty( get_thread() );
}


/* ----------------------------------------------------------------
 * Update function and event hook
 * ---------------------------------------------------------------- */

void
nest::poisson_generator::update( Time const& T, const long from, const long to )
{
  if ( order_certified_drive::handles_thread( get_thread() ) )
  {
    return;  // EXPERIMENTAL-NOT-STOCK (order-certified step): drawn by the update phase's prepass, applied by the targets
  }

  if ( P_.rate_ <= 0 )
  {
    return;
  }

  for ( long lag = from; lag < to; ++lag )
  {
    if ( not StimulationDevice::is_active( T + Time::step( lag ) ) )
    {
      continue;  // no spike at this lag
    }

    DSSpikeEvent se;
    kernel().event_delivery_manager.send( *this, se, lag );
  }
}

void
nest::poisson_generator::order_certified_draw( size_t tid,
  const Time& origin,
  long from,
  long to,
  const std::vector< size_t >& targets,
  std::vector< unsigned long >& counts,
  long slice_len )
{
  // EXPERIMENTAL-NOT-STOCK (order-certified step): the same tests and the same draws as update() and event_hook(), in the same order.
  if ( P_.rate_ <= 0 )
  {
    return;
  }
  for ( long lag = from; lag < to; ++lag )
  {
    if ( not StimulationDevice::is_active( origin + Time::step( lag ) ) )
    {
      continue;
    }
    for ( const size_t i : targets )
    {
      const long n_spikes = V_.poisson_dist_( get_vp_specific_rng( tid ) );
      counts[ i * static_cast< size_t >( slice_len ) + static_cast< size_t >( lag - from ) ] =
        n_spikes > 0 ? static_cast< unsigned long >( n_spikes ) : 0UL;
    }
  }
}

void
nest::poisson_generator::event_hook( DSSpikeEvent& e )
{
  long n_spikes = V_.poisson_dist_( get_vp_specific_rng( get_thread() ) );

  if ( n_spikes > 0 )  // we must not send events with multiplicity 0
  {
    e.set_multiplicity( n_spikes );
    e.get_receiver().handle( e );
  }
}

/* ----------------------------------------------------------------
 * Other functions
 * ---------------------------------------------------------------- */

void
nest::poisson_generator::set_data_from_stimulation_backend( std::vector< double >& input_param )
{
  Parameters_ ptmp = P_;  // temporary copy in case of errors

  // For the input backend
  if ( not input_param.empty() )
  {
    if ( input_param.size() != 1 )
    {
      throw BadParameterValue( "The size of the data for the poisson generator needs to be 1 [rate]." );
    }
    Dictionary d;
    d[ names::rate ] = input_param[ 0 ];
    ptmp.set( d, this );
  }

  // if we get here, temporary contains consistent set of properties
  P_ = ptmp;
}
