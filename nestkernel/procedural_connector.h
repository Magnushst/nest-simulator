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
 * EXPERIMENTAL, NOT STOCK NEST (regenerated connections, Stage A; this branch, EXPERIMENTAL-NOT-STOCK-procedural-
 * invariant: one realisation at any number of virtual processes, NEST_INVARIANT.md).
 *
 * Projections of NEST's fixed_total_number rule with static synapses, regenerated instead of stored. The realisation is
 * keyed by the network alone, never by a process:
 * - each source's count, Multinomial(N; 1/n) over the n sources, from one share tree keyed by the projection;
 * - a source's count split over procedural_classes classes of target (the node id modulo procedural_classes), in
 *   proportion to the classes' targets, by a binomial tree keyed by (projection, tree node, source);
 * - connection i of class c: its target uniform over the class's targets in node-id order, its weight rounded to a
 *   multiple of 2^-18 pA (so every input sum is exact in any order) and its delay, keyed by (projection, source, c, i).
 * NEST places node n on virtual process n mod P, so with P dividing procedural_classes each class lies on one process,
 * which draws only its own classes. Each process keeps, per projection, one record: its key, sources, every class's
 * size, its own classes' targets, the laws and its receptor port. Each source with connections on a process has one
 * entry in the source table and in this connector, under the synapse type "static_synapse_procedural", so NEST's target
 * tables and compressed spikes route its events as for stored connections; the entry holds the source's whole count.
 * When an event arrives, the connector splits that count over the process's classes and regenerates their connections.
 * A process with one or two classes (512 or more virtual processes) is split once when building instead: its entries
 * hold the source's count in each of its classes, so delivery walks no tree; the counts are the same either way.
 *
 * Refused with NotImplemented: reading or setting single connections, plasticity, disconnection, devices.
 */

#ifndef PROCEDURAL_CONNECTOR_H
#define PROCEDURAL_CONNECTOR_H

// C++ includes:
#include <array>
#include <cmath>
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

//! The classes of target, by node id modulo this number; a number of virtual processes must divide it.
constexpr std::uint32_t procedural_classes = 1024;
//! The share tree's process field for the one tree every process walks (no process's own).
constexpr std::uint64_t procedural_all_processes = 0xFFFFFF;
//! Weights are rounded to multiples of 1 / procedural_weight_scale pA.
constexpr double procedural_weight_scale = 262144.0;

//! One regenerated projection on one virtual process.
struct ProceduralRecord
{
  ProceduralRecord( const ProceduralWords::Philox::key_type& key_,
    const std::uint64_t vp_,
    NodeCollectionPTR sources_,
    const std::vector< std::uint64_t >& class_size,
    std::vector< std::uint32_t > own_classes,
    std::vector< std::vector< Node* > > own_targets,
    const std::uint64_t total,
    const ProceduralLaw& weight_,
    const ProceduralLaw& delay_,
    const size_t rport_ )
    : key( key_ )
    , key32( { static_cast< std::uint32_t >( key_[ 0 ] ), static_cast< std::uint32_t >( key_[ 1 ] ) } )
    , vp( vp_ )
    , sources( std::move( sources_ ) )
    , class_prefix( procedural_classes + 1, 0 )
    , own_prefix( procedural_classes + 1, 0 )
    , own_index( procedural_classes, -1 )
    , own( own_classes )
    , direct( own_classes.size() <= 2 )
    , targets( std::move( own_targets ) )
    , shares( key_, procedural_all_processes, sources->size(), total )
    , weight( weight_ )
    , delay( delay_ )
    , rport( rport_ )
  {
    if ( ( key_[ 0 ] >> 32 ) != 0 or ( key_[ 1 ] >> 32 ) != 0 or ( vp_ >> 24 ) != 0 or ( sources->size() >> 32 ) != 0 )
    {
      throw KernelException(
        "Regenerated connections need a seed and projection number below 2^32, fewer than 2^24 virtual processes "
        "and fewer than 2^32 sources per projection." );
    }
    for ( std::uint32_t c = 0; c < procedural_classes; ++c )
    {
      class_prefix[ c + 1 ] = class_prefix[ c ] + class_size[ c ];
    }
    for ( std::size_t j = 0; j < own_classes.size(); ++j )
    {
      own_index[ own_classes[ j ] ] = static_cast< std::int32_t >( j );
    }
    for ( std::uint32_t c = 0; c < procedural_classes; ++c )
    {
      own_prefix[ c + 1 ] = own_prefix[ c ] + ( own_index[ c ] >= 0 ? 1 : 0 );
    }
  }

  /**
   * Splits source s's k connections over the classes of target and calls use( class, count ) for each of this
   * process's classes with connections, in increasing order of class. The tree is the same on every process: node v
   * over the classes [lo, hi) splits at mid into 2v and 2v + 1, the left count Binomial(count, size of [lo, mid) /
   * size of [lo, hi)) keyed by (projection, v, s); a process descends only into branches holding classes of its own.
   */
  template < typename F >
  void
  class_counts( const std::uint64_t s, const std::uint64_t k, F&& use ) const
  {
    struct Branch
    {
      std::uint64_t v, count;
      std::uint32_t lo, hi;
    };
    std::array< Branch, 64 > stack;
    std::size_t top = 0;
    stack[ top++ ] = Branch { 1, k, 0, procedural_classes };
    while ( top > 0 )
    {
      const Branch b = stack[ --top ];
      if ( b.count == 0 or own_prefix[ b.hi ] == own_prefix[ b.lo ] )
      {
        continue;
      }
      if ( b.hi - b.lo == 1 )
      {
        use( b.lo, b.count );
        continue;
      }
      const std::uint32_t mid = b.lo + ( b.hi - b.lo ) / 2;
      const std::uint64_t in_left = class_prefix[ mid ] - class_prefix[ b.lo ];
      const std::uint64_t in_all = class_prefix[ b.hi ] - class_prefix[ b.lo ];
      std::uint64_t left = 0;
      if ( in_left == in_all )
      {
        left = b.count;
      }
      else if ( in_left > 0 )
      {
        ProceduralWords w( key, procedural_counter( ProceduralStream::classes, b.v, s, 0 ) );
        left = procedural_binomial( b.count, static_cast< double >( in_left ) / static_cast< double >( in_all ), w );
      }
      stack[ top++ ] = Branch { 2 * b.v + 1, b.count - left, mid, b.hi }; // the right half after the left
      stack[ top++ ] = Branch { 2 * b.v, left, b.lo, mid };
    }
  }

  //! Puts connection i of class c of source s32 into lane l: its counter for blocks 0 and 1, the key.
  void
  set_lane( ProceduralLanes& a,
    ProceduralLanes& b,
    ProceduralKeyLanes& keys,
    const std::size_t l,
    const std::uint32_t s32,
    const std::uint32_t c,
    const std::uint32_t i ) const
  {
    const ProceduralCtr32 ctr = procedural_connection_counter( s32, c, i, 0 );
    for ( std::size_t j = 0; j < 4; ++j )
    {
      a.w[ j ][ l ] = ctr[ j ];
      b.w[ j ][ l ] = ctr[ j ];
    }
    b.w[ 3 ][ l ] = 1;
    keys.w[ 0 ][ l ] = key32[ 0 ];
    keys.w[ 1 ][ l ] = key32[ 1 ];
  }

  /**
   * From lane l's two blocks, connection i of class c of source s32: its target, uniform over the class's targets in
   * node-id order, then its weight (rounded to a multiple of 2^-18 pA) and delay in ms by the laws. Delivery and the
   * delay pass both draw here.
   */
  void
  draw( const ProceduralZiggurat& z,
    const ProceduralLanes& a,
    const ProceduralLanes& b,
    const std::size_t l,
    const std::uint32_t s32,
    const std::uint32_t c,
    const std::uint32_t i,
    Node*& t,
    double& wt,
    double& d ) const
  {
    ProceduralWords32 w( key32,
      procedural_connection_counter( s32, c, i, 0 ),
      { a.w[ 0 ][ l ], a.w[ 1 ][ l ], a.w[ 2 ][ l ], a.w[ 3 ][ l ], b.w[ 0 ][ l ], b.w[ 1 ][ l ], b.w[ 2 ][ l ], b.w[ 3 ][ l ] } );
    const std::vector< Node* >& in_class = targets[ static_cast< std::size_t >( own_index[ c ] ) ];
    t = in_class[ w.below( in_class.size() ) ];
    wt = std::nearbyint( weight.draw( z, w ) * procedural_weight_scale ) / procedural_weight_scale;
    d = delay.draw( z, w );
  }

  /**
   * Regenerates source s's connections on this process (its k connections split over the classes), eight at a time
   * (the delay pass when building), calling on_class( class, count ) for each of the process's classes with connections
   * and use( target, weight, delay in ms ) for each connection; returns how many there were.
   */
  template < typename F, typename C >
  std::uint64_t
  connections( const ProceduralZiggurat& z, const std::uint64_t s, const std::uint64_t k, F&& use, C&& on_class ) const
  {
    ProceduralLanes a {};
    ProceduralLanes b {};
    ProceduralKeyLanes keys {};
    std::array< Node*, procedural_lanes > t;
    std::array< double, procedural_lanes > wt;
    std::array< double, procedural_lanes > d;
    std::array< std::pair< std::uint32_t, std::uint32_t >, procedural_lanes > ci;
    std::size_t n = 0;
    const auto s32 = static_cast< std::uint32_t >( s );
    const auto flush = [ & ]()
    {
      procedural_philox_lanes( a, keys );
      procedural_philox_lanes( b, keys );
      for ( std::size_t l = 0; l < n; ++l )
      {
        draw( z, a, b, l, s32, ci[ l ].first, ci[ l ].second, t[ l ], wt[ l ], d[ l ] );
      }
      for ( std::size_t l = 0; l < n; ++l )
      {
        use( t[ l ], wt[ l ], d[ l ] );
      }
      n = 0;
    };
    std::uint64_t here = 0;
    class_counts( s,
      k,
      [ & ]( const std::uint32_t c, const std::uint64_t kc )
      {
        on_class( c, kc );
        for ( std::uint64_t i = 0; i < kc; ++i )
        {
          set_lane( a, b, keys, n, s32, c, static_cast< std::uint32_t >( i ) );
          ci[ n ] = { c, static_cast< std::uint32_t >( i ) };
          if ( ++n == procedural_lanes )
          {
            flush();
          }
        }
        here += kc;
      } );
    if ( n > 0 )
    {
      flush();
    }
    return here;
  }

  ProceduralWords::Philox::key_type key;       //!< (seed, projection), for the trees' Philox4x64
  ProceduralKey32 key32;                       //!< the same, for the connections' Philox4x32-10
  std::uint64_t vp;                            //!< the virtual process
  NodeCollectionPTR sources;                   //!< the projection's sources, by index
  std::vector< std::uint64_t > class_prefix;   //!< the projection's targets in the classes below c, over all processes
  std::vector< std::uint32_t > own_prefix;     //!< this process's classes below c
  std::vector< std::int32_t > own_index;       //!< a class's place in targets, or -1 if not this process's
  std::vector< std::uint32_t > own;            //!< this process's classes with targets, in increasing order
  bool direct;                                 //!< one or two classes here: entries hold each class's count
  std::vector< std::vector< Node* > > targets; //!< this process's classes' targets, each in node-id order
  ProceduralShareTree shares;                  //!< each source's connections, the same tree on every process
  ProceduralLaw weight;
  ProceduralLaw delay;
  size_t rport; //!< receptor port, as the connection check gave it
};

//! The connections of one source in one record on one thread: one local connection index.
struct ProceduralEntry
{
  std::uint32_t record;
  //! the source's whole count in the projection, which delivery splits over classes; when the record is direct, its
  //! count in the process's second class (0 with one class)
  std::uint32_t count : 31;
  std::uint32_t more_targets : 1; //!< the next entry has the same source (set by the source table)
  std::uint32_t source;           //!< index among the record's sources (fewer than 2^32, checked when building)
  std::uint32_t first;            //!< when the record is direct, the source's count in the process's first class
};
static_assert( sizeof( ProceduralEntry ) == 16, "a regenerated entry stays 16 bytes" );

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

  //! Adds source's entry: its whole count, or for a direct record its counts in the process's second and first class.
  void
  add_entry( const std::uint32_t r, const std::uint64_t source, const std::uint64_t count, const std::uint64_t first )
  {
    if ( count >= ( std::uint64_t { 1 } << 31 ) or first >= ( std::uint64_t { 1 } << 31 ) )
    {
      throw KernelException( "A regenerated entry holds at most 2^31 - 1 connections per count." );
    }
    entries_.push_back( ProceduralEntry { r,
      static_cast< std::uint32_t >( count ),
      0,
      static_cast< std::uint32_t >( source ),
      static_cast< std::uint32_t >( first ) } );
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
   * The batches of eight are filled across those entries, in order, each lane with its own record's key and counter:
   * the same draws as ProceduralRecord::connections, handed to the targets in the same order, but no partly empty batch
   * per entry when entries hold few connections.
   */
  size_t
  send( const size_t, const size_t lcid, const std::vector< ConnectorModel* >&, Event& e ) override
  {
    struct Lane
    {
      const ProceduralRecord* r;
      size_t lcid;
      std::uint32_t s;
      std::uint32_t c;
      std::uint32_t i;
    };
    std::array< Lane, procedural_lanes > lane;
    ProceduralLanes a {};
    ProceduralLanes b {};
    ProceduralKeyLanes keys {};
    std::array< Node*, procedural_lanes > t;
    std::array< double, procedural_lanes > wt;
    std::array< double, procedural_lanes > d;
    std::size_t n = 0;
    const auto deliver = [ & ]()
    {
      procedural_philox_lanes( a, keys );
      procedural_philox_lanes( b, keys );
      for ( std::size_t l = 0; l < n; ++l )
      {
        lane[ l ].r->draw( zig_, a, b, l, lane[ l ].s, lane[ l ].c, lane[ l ].i, t[ l ], wt[ l ], d[ l ] );
      }
      for ( std::size_t l = 0; l < n; ++l )
      {
        e.set_port( lane[ l ].lcid );
        e.set_weight( wt[ l ] );
        e.set_delay_steps( Time::delay_ms_to_steps( d[ l ] ) );
        e.set_receiver( *t[ l ] );
        e.set_rport( lane[ l ].r->rport );
        e();
      }
      n = 0;
    };
    size_t delivered = 0;
    for ( size_t lc = lcid;; ++lc )
    {
      const ProceduralEntry& entry = entries_[ lc ];
      const ProceduralRecord& r = *records_[ entry.record ];
      const std::uint32_t s32 = entry.source;
      const auto per_class = [ & ]( const std::uint32_t c, const std::uint64_t kc )
      {
        for ( std::uint64_t i = 0; i < kc; ++i )
        {
          r.set_lane( a, b, keys, n, s32, c, static_cast< std::uint32_t >( i ) );
          lane[ n ] = { &r, lc, s32, c, static_cast< std::uint32_t >( i ) };
          if ( ++n == procedural_lanes )
          {
            deliver();
          }
        }
        delivered += kc;
      };
      if ( r.direct )
      {
        // the counts split when building, in increasing order of class as the tree gives them
        per_class( r.own[ 0 ], entry.first );
        if ( r.own.size() > 1 )
        {
          per_class( r.own[ 1 ], entry.count );
        }
      }
      else
      {
        r.class_counts( s32, entry.count, per_class );
      }
      if ( not entry.more_targets )
      {
        break;
      }
    }
    if ( n > 0 )
    {
      deliver();
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
