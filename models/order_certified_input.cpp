/*
 *  order_certified_input.cpp
 *
 *  EXPERIMENTAL, NOT STOCK NEST 3.10 (branch EXPERIMENTAL-NOT-STOCK-order-certified-step). See the header.
 */

#include "order_certified_input.h"

// C++ includes:
#include <cstdlib>
#include <iostream>
#include <mutex>

// Includes from nestkernel:
#include "kernel_manager.h"
#include "simulation_manager.h"

// Includes from models:
#include "iaf_psc_delta.h"
#include "iaf_psc_exp.h"
#include "multimeter.h"
#include "parrot_neuron.h"
#include "poisson_generator.h"
#include "spike_generator.h"
#include "spike_recorder.h"

namespace nest
{
namespace order_certified_input
{
namespace
{

std::vector< ThreadInput > tables_;
std::once_flag sized_;

//! NEST_ORDER_CERTIFIED_U2_OFF=1 keeps every neuron's input in its own buffer (the stock layout); read once.
bool
switched_off()
{
  static const bool off = [] {
    const char* v = std::getenv( "NEST_ORDER_CERTIFIED_U2_OFF" );
    return v != nullptr and v[ 0 ] == '1';
  }();
  return off;
}

bool
reporting()
{
  static const bool report = [] {
    const char* v = std::getenv( "NEST_ORDER_CERTIFIED_REPORT" );
    return v != nullptr and v[ 0 ] == '1';
  }();
  return report;
}

void
report( size_t tid, const ThreadInput& t, const char* reason )
{
  if ( reporting() )
  {
#pragma omp critical( order_certified_report )
    std::cerr << "order_certified_input: rank " << kernel().mpi_manager.get_rank() << " thread " << tid << " built "
              << t.built << " neurons " << t.n << " slots " << t.n_slots << " (" << reason << ")" << std::endl;
  }
}

//! Build thread tid's table from its iaf_psc_exp neurons' own buffers, unless a condition fails.
void
build( size_t tid, ThreadInput& t )
{
  t = ThreadInput();
  t.tried = true;
  if ( switched_off() )
  {
    report( tid, t, "switched off" );
    return;
  }
  const SparseNodeArray& nodes = kernel().node_manager.get_local_nodes( tid );
  for ( auto it = nodes.begin(); it != nodes.end(); ++it )
  {
    Node* node = it->get_node();
    if ( auto* n = dynamic_cast< iaf_psc_exp* >( node ) )
    {
      t.neurons.push_back( n );
    }
    else if ( dynamic_cast< iaf_psc_delta* >( node ) == nullptr and dynamic_cast< parrot_neuron* >( node ) == nullptr
      and dynamic_cast< poisson_generator* >( node ) == nullptr and dynamic_cast< spike_generator* >( node ) == nullptr
      and dynamic_cast< spike_recorder* >( node ) == nullptr and dynamic_cast< multimeter* >( node ) == nullptr )
    {
      // e.g. a dc_generator: its CurrentEvents would reach the neurons' own buffers, not the table
      t = ThreadInput();
      t.tried = true;
      report( tid, t, "a node that may send current events, or of another model, on the thread" );
      return;
    }
  }
  if ( t.neurons.empty() )
  {
    t = ThreadInput();
    t.tried = true;
    report( tid, t, "no iaf_psc_exp on the thread" );
    return;
  }
  t.n = t.neurons.size();
  t.n_slots = t.neurons[ 0 ]->order_certified_input_slots();
  for ( const iaf_psc_exp* n : t.neurons )
  {
    if ( n->order_certified_input_slots() != t.n_slots )
    {
      t = ThreadInput();
      t.tried = true;
      report( tid, t, "neurons with buffers of different sizes" );
      return;
    }
  }
  t.data.resize( t.n_slots * t.n );
  for ( size_t i = 0; i < t.n; ++i )
  {
    for ( size_t s = 0; s < t.n_slots; ++s )
    {
      slot( t, s, i ) = t.neurons[ i ]->order_certified_input_get( s );
    }
    t.neurons[ i ]->set_order_certified_input_index( static_cast< long >( i ) );
  }
  t.built = true;
  report( tid, t, "all conditions hold" );
}

//! Copy the table back into the neurons' own buffers, return them to the stock path, and release the table.
void
write_back( ThreadInput& t )
{
  if ( t.built )
  {
    for ( size_t i = 0; i < t.n; ++i )
    {
      for ( size_t s = 0; s < t.n_slots; ++s )
      {
        t.neurons[ i ]->order_certified_input_set( s, slot( t, s, i ) );
      }
      t.neurons[ i ]->set_order_certified_input_index( -1 );
    }
  }
  t = ThreadInput();
}

void
post_run_hook()
{
#pragma omp parallel for schedule( static )
  for ( size_t tid = 0; tid < tables_.size(); ++tid )
  {
    write_back( tables_[ tid ] );
  }
}

}  // namespace

void
register_hooks()
{
  std::call_once( sized_, [] { tables_.resize( kernel().vp_manager.get_num_threads() ); } );
  order_certified_post_run.store( &post_run_hook, std::memory_order_release );
}

ThreadInput&
table( size_t tid )
{
  return tables_[ tid ];
}

void
prepass( size_t tid )
{
  if ( tid < tables_.size() and not tables_[ tid ].tried )
  {
    build( tid, tables_[ tid ] );
  }
}

}  // namespace order_certified_input
}  // namespace nest
