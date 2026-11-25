/*
 *  order_certified_drive.cpp
 *
 *  EXPERIMENTAL, NOT STOCK NEST 3.10 (branch EXPERIMENTAL-NOT-STOCK-order-certified-step). See the header.
 */

#include "order_certified_drive.h"

// C++ includes:
#include <mutex>
#include <unordered_map>

// Includes from nestkernel:
#include "connector_base.h"
#include "kernel_manager.h"
#include "simulation_manager.h"
#include "target_identifier.h"

// Includes from models:
#include "dc_generator.h"
#include "iaf_psc_delta.h"
#include "iaf_psc_exp.h"
#include "multimeter.h"
#include "parrot_neuron.h"
#include "poisson_generator.h"
#include "spike_generator.h"
#include "spike_recorder.h"
#include "static_synapse.h"

namespace nest
{
namespace order_certified_drive
{
namespace
{

std::vector< ThreadTable > tables_;
std::once_flag sized_;

//! Append the targets of one connector of a generator; false if a connection is not eligible.
template < typename TargetIdentifierT >
bool
append_targets( const ConnectorBase* cb,
  size_t tid,
  std::vector< size_t >& gen_targets,
  std::vector< Node* >& targets,
  std::vector< double >& weight,
  std::vector< long >& delay,
  std::unordered_map< Node*, size_t >& index )
{
  const auto* c = dynamic_cast< const Connector< static_synapse< TargetIdentifierT > >* >( cb );
  if ( c == nullptr )
  {
    return false;
  }
  const auto& conns = c->get_connections_order_certified();
  for ( size_t lcid = 0; lcid < conns.size(); ++lcid )  // the order of Connector::send_to_all
  {
    const auto& conn = conns[ lcid ];
    Node* target = conn.get_target( tid );
    if ( dynamic_cast< iaf_psc_exp* >( target ) == nullptr and dynamic_cast< iaf_psc_delta* >( target ) == nullptr )
    {
      return false;  // only these two neuron models apply the drive in their update
    }
    if ( index.count( target ) != 0 )
    {
      return false;  // a second drive connection to one target: keep the stock path
    }
    if ( target->is_frozen() )
    {
      return false;  // stock adds to a frozen target's buffer, but its update() never runs here
    }
    const size_t i = targets.size();
    index[ target ] = i;
    targets.push_back( target );
    weight.push_back( conn.get_weight_order_certified() );
    delay.push_back( conn.get_delay_steps() );
    gen_targets.push_back( i );
  }
  return true;
}

//! Build thread tid's table; leave it disabled, with every neuron on the stock path, unless all conditions hold.
void
build( size_t tid, ThreadTable& tt )
{
  tt = ThreadTable();
  tt.dirty = false;

  std::vector< Node* > targets;
  std::unordered_map< Node*, size_t > index;
  bool ok = true;

  // 1. Every node on the thread is of a known model, and none but the generators draws in the update phase.
  const SparseNodeArray& nodes = kernel().node_manager.get_local_nodes( tid );
  for ( auto it = nodes.begin(); it != nodes.end(); ++it )
  {
    Node* node = it->get_node();
    if ( auto* n = dynamic_cast< iaf_psc_exp* >( node ) )
    {
      n->set_order_certified_index( -1 );
      ok = ok and n->order_certified_deterministic();  // no escape noise, so no draw in its update
    }
    else if ( auto* n = dynamic_cast< iaf_psc_delta* >( node ) )
    {
      n->set_order_certified_index( -1 );
    }
    else if ( auto* g = dynamic_cast< poisson_generator* >( node ) )
    {
      tt.generators.push_back( g );
    }
    else if ( dynamic_cast< spike_generator* >( node ) != nullptr )
    {
      ok = false;  // it could add to a neuron's slot in the update phase, in between the drive
    }
    else if ( dynamic_cast< parrot_neuron* >( node ) == nullptr and dynamic_cast< dc_generator* >( node ) == nullptr
      and dynamic_cast< spike_recorder* >( node ) == nullptr and dynamic_cast< multimeter* >( node ) == nullptr )
    {
      ok = false;  // an unknown model might draw from the thread's stream in the update phase
    }
  }

  // 2. Every generator's connections are static synapses to iaf_psc_exp or iaf_psc_delta, one per target.
  for ( size_t g = 0; ok and g < tt.generators.size(); ++g )
  {
    tt.gen_targets.emplace_back();
    const auto& connectors =
      kernel().connection_manager.get_connectors_from_device_order_certified( tid, tt.generators[ g ]->get_local_device_id() );
    for ( const ConnectorBase* cb : connectors )  // by synapse type, as send_from_device() walks them
    {
      if ( cb == nullptr )
      {
        continue;
      }
      if ( not append_targets< TargetIdentifierPtrRport >(
             cb, tid, tt.gen_targets.back(), targets, tt.weight, tt.delay, index )
        and not append_targets< TargetIdentifierIndex >(
          cb, tid, tt.gen_targets.back(), targets, tt.weight, tt.delay, index ) )
      {
        ok = false;
        break;
      }
    }
  }

  if ( not ok or tt.generators.empty() )
  {
    tt = ThreadTable();
    tt.dirty = false;
    return;
  }

  for ( size_t i = 0; i < targets.size(); ++i )
  {
    if ( auto* n = dynamic_cast< iaf_psc_exp* >( targets[ i ] ) )
    {
      n->set_order_certified_index( static_cast< long >( i ) );
    }
    else
    {
      dynamic_cast< iaf_psc_delta* >( targets[ i ] )->set_order_certified_index( static_cast< long >( i ) );
    }
  }
  tt.enabled = true;
}

void
prepass_hook( size_t tid, const Time& origin, long from, long to )
{
  prepass( tid, origin, from, to );
}

}  // namespace

void
register_hook()
{
  std::call_once( sized_, [] { tables_.resize( kernel().vp_manager.get_num_threads() ); } );
  order_certified_update_prepass.store( &prepass_hook, std::memory_order_release );
}

void
mark_dirty( size_t tid )
{
  if ( tid < tables_.size() )
  {
    tables_[ tid ].dirty = true;
  }
}

ThreadTable&
table( size_t tid )
{
  return tables_[ tid ];
}

bool
handles_thread( size_t tid )
{
  return tid < tables_.size() and not tables_[ tid ].dirty and tables_[ tid ].enabled;
}

void
prepass( size_t tid, const Time& origin, long from, long to )
{
  if ( tid >= tables_.size() )
  {
    return;  // more threads than when the tables were sized: this thread keeps the stock path
  }
  ThreadTable& tt = tables_[ tid ];
  if ( tt.dirty )
  {
    build( tid, tt );
  }
  if ( not tt.enabled )
  {
    return;
  }
  tt.slice_len = to - from;
  tt.counts.assign( tt.weight.size() * static_cast< size_t >( tt.slice_len ), 0 );
  for ( size_t g = 0; g < tt.generators.size(); ++g )
  {
    poisson_generator* gen = tt.generators[ g ];
    if ( gen->is_frozen() )
    {
      continue;  // the node loop does not update frozen nodes
    }
    gen->order_certified_draw( tid, origin, from, to, tt.gen_targets[ g ], tt.counts, tt.slice_len );
  }
}

}  // namespace order_certified_drive
}  // namespace nest
