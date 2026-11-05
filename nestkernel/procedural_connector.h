/*
 *  procedural_connector.h
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
 * Projections of NEST's fixed_total_number rule with static synapses, regenerated instead of stored. Each virtual
 * process keeps, per projection, one record: its key, sources, local targets, its connections' count per source (a
 * share tree, walked once when building), its weight and delay laws and its receptor port. Each source with connections
 * on a process has one entry in the source table and in this connector, under the synapse type
 * "static_synapse_procedural", so NEST's target tables and compressed spikes route its events as for stored
 * connections; the entry holds that source's count. When an event arrives, the connector regenerates that many
 * connections from the record's key and hands each event to its target as static_synapse does. The same key gives the
 * same connections at every delivery.
 *
 * Refused with NotImplemented: reading or setting single connections, plasticity, disconnection, devices.
 */

#ifndef PROCEDURAL_CONNECTOR_H
#define PROCEDURAL_CONNECTOR_H

// C++ includes:
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Includes from libnestutil:
#include "block_vector.h"
#include "sort.h"

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connector_base.h"
#include "connector_model.h"
#include "event.h"
#include "exceptions.h"
#include "nest_names.h"
#include "nest_time.h"
#include "node.h"
#include "node_collection.h"
#include "procedural_random.h"
#include "procedural_share.h"

namespace nest
{

//! One regenerated projection on one virtual process.
struct ProceduralRecord
{
  ProceduralRecord( const ProceduralWords::Philox::key_type& key_,
    const std::uint64_t vp_,
    NodeCollectionPTR sources_,
    std::vector< Node* > targets_,
    const std::uint64_t total,
    const ProceduralLaw& weight_,
    const ProceduralLaw& delay_,
    const size_t rport_ )
    : key( key_ )
    , vp( vp_ )
    , sources( std::move( sources_ ) )
    , targets( std::move( targets_ ) )
    , shares( key_, vp_, sources->size(), total )
    , weight( weight_ )
    , delay( delay_ )
    , rport( rport_ )
  {
  }

  /**
   * Regenerates connection i of source s on this process and calls use( target index, weight, delay in ms ). The
   * target is uniform over the process's targets, as in NEST's fixed_total_number; weight and delay follow the laws.
   */
  template < typename F >
  void
  connection( const ProceduralZiggurat& z, const std::uint64_t s, const std::uint64_t i, F&& use ) const
  {
    ProceduralWords w( key, procedural_counter( ProceduralStream::connection, i, s, vp ) );
    const std::uint64_t t = w.below( targets.size() );
    const double wt = weight.draw( z, w );
    const double d = delay.draw( z, w );
    use( t, wt, d );
  }

  ProceduralWords::Philox::key_type key; //!< (seed, projection)
  std::uint64_t vp;                      //!< the virtual process
  NodeCollectionPTR sources;             //!< the projection's sources, by index
  std::vector< Node* > targets;          //!< its targets on this process, in NEST's order
  ProceduralShareTree shares;            //!< its connections per source on this process
  ProceduralLaw weight;
  ProceduralLaw delay;
  size_t rport; //!< receptor port, as the connection check gave it
};

//! The connections of one source in one record on one thread: one local connection index.
struct ProceduralEntry
{
  std::uint32_t record;
  std::uint32_t count : 31;       //!< how many, from the share tree when building, so that delivery walks no tree
  std::uint32_t more_targets : 1; //!< the next entry has the same source (set by the source table)
  std::uint64_t source;           //!< index among the record's sources
};

class ProceduralConnector : public ConnectorBase
{
public:
  explicit ProceduralConnector( const synindex syn_id )
    : syn_id_( syn_id )
  {
  }

  std::uint32_t
  add_record( std::unique_ptr< ProceduralRecord > r )
  {
    records_.push_back( std::move( r ) );
    return static_cast< std::uint32_t >( records_.size() - 1 );
  }

  const ProceduralRecord&
  record( const std::uint32_t r ) const
  {
    return *records_[ r ];
  }

  void
  add_entry( const std::uint32_t r, const std::uint64_t source, const std::uint64_t count )
  {
    if ( count >= ( std::uint64_t { 1 } << 31 ) )
    {
      throw KernelException( "A regenerated entry holds at most 2^31 - 1 connections." );
    }
    entries_.push_back( ProceduralEntry { r, static_cast< std::uint32_t >( count ), 0, source } );
  }

  synindex
  get_syn_id() const override
  {
    return syn_id_;
  }

  size_t
  size() const override
  {
    return entries_.size();
  }

  /**
   * Delivers the event to every regenerated connection of the entries from lcid on that share its source, as
   * static_synapse::send does for each, and returns how many connections it reached (the cycle log counts them).
   */
  size_t
  send( const size_t, const size_t lcid, const std::vector< ConnectorModel* >&, Event& e ) override
  {
    size_t delivered = 0;
    for ( size_t l = lcid;; ++l )
    {
      const ProceduralEntry& entry = entries_[ l ];
      const ProceduralRecord& r = *records_[ entry.record ];
      const std::uint64_t k = entry.count;
      e.set_port( l );
      for ( std::uint64_t i = 0; i < k; ++i )
      {
        r.connection( zig_,
          entry.source,
          i,
          [ & ]( const std::uint64_t t, const double weight, const double delay )
          {
            e.set_weight( weight );
            e.set_delay_steps( Time::delay_ms_to_steps( delay ) );
            e.set_receiver( *r.targets[ t ] );
            e.set_rport( r.rport );
            e();
          } );
      }
      delivered += k;
      if ( not entry.more_targets )
      {
        break;
      }
    }
    return delivered;
  }

  void
  sort_connections( BlockVector< Source >& sources ) override
  {
    nest::sort( sources, entries_ );
  }

  void
  set_source_has_more_targets( const size_t lcid, const bool has_more_targets ) override
  {
    entries_[ lcid ].more_targets = has_more_targets;
  }

  //! Static synapses do not change.
  void
  send_weight_event( const size_t, const unsigned int, Event&, const CommonSynapseProperties& ) override
  {
  }

  //! Static synapses do not change.
  void
  trigger_update_weight( const long,
    const size_t,
    const std::vector< spikecounter >&,
    const double,
    const std::vector< ConnectorModel* >& ) override
  {
  }

  void
  get_synapse_status( const size_t, const size_t, Dictionary& ) const override
  {
    refuse_( "Reading single connections" );
  }

  void
  set_synapse_status( const size_t, const Dictionary&, ConnectorModel& ) override
  {
    refuse_( "Setting single connections" );
  }

  void
  get_connection( const size_t,
    const size_t,
    const size_t,
    const size_t,
    const long,
    std::deque< ConnectionID >& ) const override
  {
    refuse_( "GetConnections" );
  }

  void
  get_connection_with_specified_targets( const size_t,
    const std::vector< size_t >&,
    const size_t,
    const size_t,
    const long,
    std::deque< ConnectionID >& ) const override
  {
    refuse_( "GetConnections" );
  }

  void
  get_all_connections( const size_t, const size_t, const size_t, const long, std::deque< ConnectionID >& ) const override
  {
    refuse_( "GetConnections" );
  }

  void
  get_source_lcids( const size_t, const size_t, std::vector< size_t >& ) const override
  {
    refuse_( "Looking up connections by target" );
  }

  void
  get_target_node_ids( const size_t, const size_t, const std::string&, std::vector< size_t >& ) const override
  {
    refuse_( "Structural plasticity" );
  }

  size_t
  get_target_node_id( const size_t, const unsigned int ) const override
  {
    refuse_( "Looking up a connection's target" );
  }

  void
  send_to_all( const size_t, const std::vector< ConnectorModel* >&, Event& ) override
  {
    refuse_( "Connections from devices" );
  }

  size_t
  find_first_target( const size_t, const size_t, const size_t ) const override
  {
    refuse_( "Disconnect" );
  }

  size_t
  find_enabled_connection( const size_t, const size_t, const size_t, const size_t, const SourceTable& ) const override
  {
    refuse_( "Disconnect" );
  }

  void
  disable_connection( const size_t ) override
  {
    refuse_( "Disconnect" );
  }

  void
  remove_disabled_connections( const size_t ) override
  {
    refuse_( "Disconnect" );
  }

private:
  [[noreturn]] static void
  refuse_( const std::string& what )
  {
    throw NotImplemented( what + " is not available for regenerated connections (static_synapse_procedural)." );
  }

  synindex syn_id_;
  std::vector< std::unique_ptr< ProceduralRecord > > records_;
  BlockVector< ProceduralEntry > entries_;
  ProceduralZiggurat zig_;
};

//! The synapse type of regenerated connections. Its connections are made only by fixed_total_number with procedural.
class ProceduralConnectorModel : public ConnectorModel
{
public:
  explicit ProceduralConnectorModel( const std::string name )
    : ConnectorModel( name, ConnectionModelProperties::IS_PRIMARY | ConnectionModelProperties::HAS_DELAY )
  {
  }

  void
  add_connection( Node&, Node&, std::vector< ConnectorBase* >&, const synindex, const Dictionary&, const double, const double )
    override
  {
    throw NotImplemented(
      "static_synapse_procedural connections are made only by fixed_total_number with \"procedural\": True." );
  }

  ConnectorModel*
  clone( std::string, synindex ) const override
  {
    throw NotImplemented( "CopyModel is not available for static_synapse_procedural." );
  }

  void
  calibrate( const TimeConverter& ) override
  {
  }

  void
  get_status( Dictionary& d ) const override
  {
    d[ names::synapse_model ] = name_;
    d[ names::synapse_modelid ] = static_cast< long >( syn_id_ );
    d[ names::requires_symmetric ] = false;
    d[ names::has_delay ] = true;
  }

  void
  set_status( const Dictionary& ) override
  {
  }

  const CommonSynapseProperties&
  get_common_properties() const override
  {
    return cp_;
  }

  void
  check_synapse_params( const Dictionary& ) const override
  {
  }

  std::unique_ptr< SecondaryEvent >
  get_secondary_event() override
  {
    return nullptr;
  }

  size_t
  get_syn_id() const override
  {
    return syn_id_;
  }

  void
  set_syn_id( const synindex syn_id ) override
  {
    syn_id_ = syn_id;
  }

private:
  CommonSynapseProperties cp_;
  synindex syn_id_ = invalid_synindex;
};

} // namespace nest

#endif /* #ifndef PROCEDURAL_CONNECTOR_H */
