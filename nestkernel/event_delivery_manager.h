/*
 *  event_delivery_manager.h
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

#ifndef EVENT_DELIVERY_MANAGER_H
#define EVENT_DELIVERY_MANAGER_H

// C++ includes:
#include <cassert>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

// Includes from libnestutil:
#include "manager_interface.h"
#include "stopwatch.h"

// Includes from nestkernel:
#include "buffer_resize_log.h"
#include "event.h"
#include "mpi_manager.h"  // OffGridSpike
#include "nest_time.h"
#include "nest_types.h"
#include "node.h"
#include "per_thread_bool_indicator.h"
#include "secondary_event.h"
#include "spike_data.h"
#include "target_table.h"
#include "vp_manager.h"


namespace nest
{
typedef MPIManager::OffGridSpike OffGridSpike;

class TargetData;
class SendBufferPosition;
class TargetSendBufferPosition;


class EventDeliveryManager : public ManagerInterface
{
public:
  EventDeliveryManager();
  ~EventDeliveryManager() override;

  void initialize( const bool ) override;
  void finalize( const bool ) override;
  void set_status( const Dictionary& ) override;
  void get_status( Dictionary& ) override;

  /**
   * Standard routine for sending events.
   *
   * This method decides if the event has to be delivered locally or globally. It exists
   * to keep a clean and unitary interface for the event sending mechanism.
   * @note Only specializations of SpikeEvent send remotely. A specialization for
   *       DSSpikeEvent exists to avoid that these events are sent to remote processes.
   * \see send_local()
   */
  template < class EventT >
  void send( Node& source, EventT& e, const long lag = 0 );

  /**
   * Send a secondary event remote.
   */
  void send_secondary( Node& source, SecondaryEvent& e );

  /**
   * Send event e to all targets of node source on thread t
   */
  void send_local( size_t t, Node& source, Event& e );

  /**
   * Add node ID of event sender to the spike_register.
   *
   * An event sent through this method will remain in the queue until
   * the network time has advanced by min_delay_ steps. After this period
   * the buffers are collocated and sent to the partner machines.

   * Old documentation from network.h:
   * Place an event in the global event queue.
   * Add event to the queue to be delivered
   * when it is due.
   * At the delivery time, the target list of the sender is iterated
   * and the event is delivered to all targets.
   * The event is guaranteed to arrive at the receiver when all
   * elements are updated and the system is
   * in a synchronised (single threaded) state.
   * @see send_to_targets()
   */
  void send_remote( size_t tid, SpikeEvent&, const long lag = 0 );

  /**
   * EXPERIMENTAL, NOT STOCK: send a spike decided one step ahead by the fused sweep. As send< SpikeEvent >, but the
   * event carries the stamp given, and lag 0 in the slice that will deliver it.
   */
  void send_decided( Node& source, SpikeEvent& e, const Time& stamp );

  //! EXPERIMENTAL: whether the fused sweep is on (environment variable NEST_ASYNC set to fused or async).
  bool
  fused() const
  {
    return fused_;
  }

  //! EXPERIMENTAL: keep the spikes just received (stamped `stamp`, in steps) until the slice that delivers them.
  void hold_received_spikes( const long stamp );

  //! EXPERIMENTAL: make the held spikes stamped `stamp` the ones the coming slice delivers.
  void prepare_fused_delivery( const long stamp );

  //! EXPERIMENTAL: whether the asynchronous exchange is on (NEST_ASYNC=async).
  bool
  async_exchange() const
  {
    return async_;
  }

  /**
   * EXPERIMENTAL: post this slice's new spikes (stamped `stamp`) to every rank, point to point, without waiting.
   * Master thread only.
   */
  void async_send( const long stamp );

  //! EXPERIMENTAL: receive and hold every stamp in flight, and complete all sends (end of a Run call). Master only.
  void async_settle();

  //! EXPERIMENTAL: whether the lockstep arm is on (NEST_ASYNC=lockstep_p2p).
  bool
  lockstep_p2p() const
  {
    return lockstep_p2p_;
  }

  //! EXPERIMENTAL: the lockstep arm's exchange of this slice's spikes, completed before it returns. Master only.
  void lockstep_p2p_exchange( const long stamp );

  //! EXPERIMENTAL: whether the plain asynchronous exchange is on (NEST_ASYNC=async_plain).
  bool
  plain_async() const
  {
    return plain_async_;
  }

  //! EXPERIMENTAL: the steps by which the stamps of held spikes are shifted at delivery (plain: minus the slice).
  void
  set_delivery_offset( const long offset )
  {
    fused_deliver_offset_ = offset;
  }

  //! EXPERIMENTAL: size the per-thread stores of the plain exchange; outside parallel regions only.
  void prepare_plain_async( const size_t num_threads );

  //! EXPERIMENTAL: whether a device's input to a neuron on thread tid is kept for the next slice.
  bool
  deferring_device_input( const size_t tid ) const
  {
    return plain_async_ and not plain_break_no_defer_ and deferring_[ tid ] != 0;
  }

  //! EXPERIMENTAL: thread tid starts or ends its update phase (the plain exchange keeps device inputs in it).
  void
  set_deferring( const size_t tid, const bool on )
  {
    deferring_[ tid ] = on ? 1 : 0;
  }

  //! EXPERIMENTAL: keep a device's input, rel steps after the current slice origin, for the next slice.
  void defer_device_input( const size_t tid, Node* target, const long rel, const int channel, const double value );

  //! EXPERIMENTAL: add the inputs thread tid kept, in the order they came.
  void apply_deferred_inputs( const size_t tid );

  //! EXPERIMENTAL: the spike count and digest of this rank, summed over threads (check 5)
  void experimental_spike_digest( std::uint64_t& count, std::uint64_t& digest ) const;

  //! EXPERIMENTAL: the master's cumulative exchange times (s) and what it posted, for the cycle log.
  void
  get_exchange_totals( double& post, double& wait, double& receive, double& complete, size_t& messages,
    size_t& bytes ) const
  {
    post = async_time_post_;
    wait = async_time_wait_;
    receive = async_time_receive_;
    complete = async_time_complete_;
    messages = async_messages_sent_;
    bytes = async_bytes_sent_;
  }

  /**
   * Add node ID of event sender to the spike_register.
   *
   * Store event offset with node ID.
   * An event sent through this method will remain in the queue until
   * the network time has advanced by min_delay_ steps. After this period
   * the buffers are collocated and sent to the partner machines.

   * Old documentation from network.h:
   * Place an event in the global event queue.
   * Add event to the queue to be delivered
   * when it is due.
   * At the delivery time, the target list of the sender is iterated
   * and the event is delivered to all targets.
   * The event is guaranteed to arrive at the receiver when all
   * elements are updated and the system is
   * in a synchronised (single threaded) state.
   * @see send_to_targets()
   */
  void send_off_grid_remote( size_t tid, SpikeEvent& e, const long lag = 0 );

  /**
   * Send event e directly to its target node.
   *
   * This should be
   * used only where necessary, e.g. if a node wants to reply
   * to a *RequestEvent immediately.
   */
  void send_to_node( Event& e );

  /**
   * return current communication style.
   * A result of true means off_grid, false means on_grid communication.
   */
  bool get_off_grid_communication() const;

  /**
   * set communication style to off_grid (true) or on_grid
   */
  void set_off_grid_communication( bool off_grid_spiking );

  /**
   * Return 0 for even, 1 for odd time slices.
   *
   * This is useful for buffers that need to be written alternatingly
   * by time slice. The value is given by get_slice_() % 2.
   * @see read_toggle
   */
  size_t write_toggle() const;

  /**
   * Return 1 - write_toggle().
   *
   * This is useful for buffers that need to be read alternatingly
   * by slice. The value is given by 1-write_toggle().
   * @see write_toggle
   */
  size_t read_toggle() const;

  /**
   * Return (T+d) mod max_delay.
   */
  long get_modulo( long d );


  /**
   * Index to slice-based buffer.
   * Return ((T+d)/min_delay) % ceil(max_delay/min_delay).
   */
  long get_slice_modulo( long d );

  /**
   * Resize spike_register and comm_buffer to correct dimensions.
   *
   * Resizes also off_grid_*_buffer_.
   * This is done by simulate() when called for the first time.
   * The spike buffers cannot be reconfigured later, whence neither
   * the number of local threads or the min_delay can change after
   * simulate() has been called. ConnectorModel::check_delay() and
   * Network::set_status() ensure this.
   */
  void configure_spike_data_buffers();

  void configure_spike_register();

  void resize_send_recv_buffers_target_data();

  void configure_secondary_buffers();

  /**
   * Collocates spikes from register to MPI buffers, communicates via
   * MPI and delivers events to targets.
   */
  void gather_spike_data();

  /**
   * Collocates presynaptic connection information, communicates via
   * MPI and creates presynaptic connection infrastructure.
   */
  void gather_target_data( const size_t tid );

  void gather_target_data_compressed( const size_t tid );


  /**
   * Delivers events to targets.
   */
  void deliver_events( const size_t tid );

  /**
   * Collocates presynaptic connection information for secondary events (MPI
   * buffer offsets), communicates via MPI and create presynaptic connection
   * infrastructure for secondary events.
   */
  void gather_secondary_target_data();

  void write_done_marker_secondary_events_( const bool done );

  void gather_secondary_events( const bool done );

  bool deliver_secondary_events( const size_t tid, const bool called_from_wfr_update );

  /**
   * Update modulo table based on current time settings.
   *
   * This function is called after all nodes have been updated.
   * We can compute the value of (T+d) mod max_delay without explicit
   * reference to the network clock, because compute_moduli_ is
   * called whenever the network clock advances.
   * The various modulos for all available delays are stored in
   * a lookup-table and this table is rotated once per time slice.
   *
   * Update table of fixed modulos, including slice-based.
   */
  void update_moduli();

  /**
   * Initialize modulo table.
   *
   * TODO: can probably be private
   */
  void init_moduli();

  /**
   * Set local spike counter to zero.
   */
  virtual void reset_counters();

  /**
   * Set time measurements for internal profiling to zero (reg. prep.)
   */
  virtual void reset_timers_for_preparation();

  /**
   * Set time measurements for internal profiling to zero (reg. sim. dyn.)
   */
  virtual void reset_timers_for_dynamics();

#ifdef CYCLE_TIMERS
  /**
   * Returns the elapsed time of sw_communicate_spike_data_
   * Used in SimulationManager::update_() for cycle timers
   */
  Stopwatch< StopwatchGranularity::Detailed, StopwatchParallelism::MasterOnly > get_sw_communicate_spike_data() const;

  /**
   * Returns the elapsed time of sw_collocate_spike_data_
   * Used in SimulationManager::update_() for cycle timers
   */
  Stopwatch< StopwatchGranularity::Detailed, StopwatchParallelism::MasterOnly > get_sw_collocate_spike_data() const;

  /**
   * Accumulates the spike counts stored in local_spike_counter_
   * Used in SimulationManager::update_() for cycle timers
   */
  size_t get_local_spike_counter() const;

  /**
   * Sum over threads of the connections visited by deliver_events(), cumulative
   * Used in SimulationManager::update_() for cycle timers
   */
  size_t get_local_deliver_connections() const;

  /**
   * Spike buffer section per rank, in entries, in the first and in the last all-to-all of the latest
   * gather_spike_data(); they differ only when the buffer had to grow
   * Used in SimulationManager::update_() for cycle timers
   */
  size_t get_spike_buffer_first() const;
  size_t get_spike_buffer_last() const;
#endif

private:
  template < typename SpikeDataT >
  void gather_spike_data_( std::vector< SpikeDataT >& send_buffer, std::vector< SpikeDataT >& recv_buffer );

  void resize_send_recv_buffers_spike_data_();

  /**
   * Moves spikes from on grid and off grid spike registers to correct
   * locations in MPI buffers.
   *
   * Accumulates count of spikes to be sent to any rank in num_spikes_per_rank.
   * Passed as argument, so that values accumulate as we call once for plain and once for offgrid spikes.
   */
  template < typename SpikeDataWithRankT, typename SpikeDataT >
  void collocate_spike_data_buffers_( SendBufferPosition& send_buffer_position,
    std::vector< std::vector< SpikeDataWithRankT >* >& spike_register,
    std::vector< SpikeDataT >& send_buffer,
    std::vector< size_t >& num_spikes_per_rank );

  /**
   * Set end marker for per-rank-chunks signalling completion and providing shrink/grow information.
   */
  template < typename SpikeDataT >
  void set_end_marker_( const SendBufferPosition& send_buffer_position,
    std::vector< SpikeDataT >& send_buffer,
    size_t local_max_spikes_per_rank );

  /**
   * Resets marker in MPI buffer that signals end of communication
   * across MPI ranks.
   */
  template < typename SpikeDataT >
  void reset_complete_marker_spike_data_( const SendBufferPosition& send_buffer_position,
    std::vector< SpikeDataT >& send_buffer ) const;

  /**
   * Get required buffer size.
   *
   * @returns maximum of required buffer sizes communicated by all ranks
   */
  template < typename SpikeDataT >
  size_t get_global_max_spikes_per_rank_( const SendBufferPosition& send_buffer_position,
    std::vector< SpikeDataT >& recv_buffer ) const;


  /**
   * Reads spikes from MPI buffers and delivers them to ringbuffer of
   * nodes.
   */
  template < typename SpikeDataT >
  void deliver_events_( const size_t tid,
    const std::vector< SpikeDataT >& recv_buffer,
    const size_t chunk_size_override = 0,
    const long stamp_offset = 0 );

  /**
   * Deletes all spikes from spike registers and resets spike
   * counters.
   */
  void reset_spike_register_( const size_t tid );

  /**
   * Returns true if spike has been moved to MPI buffer, such that it
   * can be removed by clean_spike_register. Required static function
   * by std::remove_if.
   */
  static bool is_marked_for_removal_( const Target& target );

  /**
   * Fills MPI buffer for communication of connection information from
   * presynaptic to postsynaptic side. Builds TargetData objects from
   * SourceTable and connections information.
   */
  bool collocate_target_data_buffers_( const size_t tid,
    const AssignedRanks& assigned_ranks,
    TargetSendBufferPosition& send_buffer_position );

  bool collocate_target_data_buffers_compressed_( const size_t tid,
    const AssignedRanks& assigned_ranks,
    TargetSendBufferPosition& send_buffer_position );

  /**
   * Sets marker in MPI buffer that signals end of communication
   * across MPI ranks.
   */
  void set_complete_marker_target_data_( const AssignedRanks& assigned_ranks,
    const TargetSendBufferPosition& send_buffer_position );

  /**
   * Reads TargetData objects from MPI buffers and creates Target
   * objects on TargetTable (presynaptic part of connection
   * infrastructure).
   */
  bool distribute_target_data_buffers_( const size_t tid );

  /**
   * Sends event e to all targets of node source. Delivers events from
   * devices directly to targets.
   */
  template < class EventT >
  void send_local_( Node& source, EventT& e, const long lag );
  void send_local_( Node& source, SecondaryEvent& e, const long lag );

  //--------------------------------------------------//

  bool off_grid_spiking_;  //!< indicates whether spikes are not constrained to
                           //!< the grid

  /**
   * Table of pre-computed modulos.
   *
   * This table is used to map time steps, given as offset from now,
   * to ring-buffer bins.  There are min_delay+max_delay bins in a ring buffer,
   * and the moduli_ array is rotated by min_delay elements after
   * each slice is completed.
   * @see RingBuffer
   */
  std::vector< long > moduli_;

  /**
   * Table of pre-computed slice-based modulos.
   *
   * This table is used to map time steps, give as offset from now,
   * to slice-based ring-buffer bins.  There are ceil(max_delay/min_delay)
   * bins in a slice-based ring buffer, one per slice within max_delay.
   * Since max_delay may not be a multiple of min_delay, we cannot simply
   * rotate the table content after each slice, but have to recompute
   * the table anew.
   * @see SliceRingBuffer
   */
  std::vector< long > slice_moduli_;

  /**
   * Register of emitted spikes.
   *
   * All spikes to be delivered non-locally are first written to this register by the thread generating the spike.
   * They are later transferred to communication buffers and exchanged globally.
   *
   * The outer dimension represents the thread generating the spikes, the second dimension the individual spikes.
   *
   * @note We store here pointers to the vectors for the individual threads so that those vectors, including their
   * administrative metadata will be stored in thread-local memory.
   */
  std::vector< std::vector< SpikeDataWithRank >* > emitted_spikes_register_;

  /**
   * Register of emitted off-grid spikes.
   *
   * All off-grid spikes to be delivered non-locally are first written to this register by the thread generating the
   * spike. They are later transferred to communication buffers and exchanged globally.
   *
   * The outer dimension represents the thread generating the spikes, the second dimension the individual spikes.
   *
   * @note We store here pointers to the vectors for the individual threads so that those vectors, including their
   * administrative metadata will be stored in thread-local memory.
   */
  std::vector< std::vector< OffGridSpikeDataWithRank >* > off_grid_emitted_spikes_register_;

  /**
   * Buffer to collect the secondary events after serialization.
   */
  std::vector< unsigned int > send_buffer_secondary_events_;
  std::vector< unsigned int > recv_buffer_secondary_events_;

  /**
   * Number of generated spike events (both off- and on-grid) during the last call to simulate.
   */
  std::vector< unsigned long > local_spike_counter_;

  //! Connections visited by deliver_events() on each thread, cumulative; counted only with cycle timers.
  std::vector< size_t > local_deliver_connections_;

  //! Spike buffer section per rank in the first and last all-to-all of the latest gather; cycle timers only.
  size_t spike_buffer_first_;
  size_t spike_buffer_last_;

  std::vector< SpikeData > send_buffer_spike_data_;
  std::vector< SpikeData > recv_buffer_spike_data_;

  // EXPERIMENTAL, NOT STOCK: the fused half-step sweep
  bool fused_ = false;
  //! received spike buffers held until their delivery, by stamp in steps, with the chunk size they were sent with
  std::map< long, std::pair< std::vector< SpikeData >, size_t > > fused_held_;
  //! the buffer the coming slice delivers, its chunk size, and whether there is one
  std::vector< SpikeData > fused_deliver_buffer_;
  size_t fused_deliver_chunk_ = 0;
  bool fused_deliver_ready_ = false;
  //! EXPERIMENTAL: the steps by which held spikes' stamps are shifted at delivery
  long fused_deliver_offset_ = 0;
  //! EXPERIMENTAL: the plain asynchronous exchange, its negative control, and its kept device inputs
  bool plain_async_ = false;
  bool plain_break_no_defer_ = false;
  std::vector< char > deferring_;
  struct DeferredInput
  {
    Node* target;
    long abs_step;
    int channel;
    double value;
  };
  std::vector< std::vector< DeferredInput > > deferred_inputs_;
  //! EXPERIMENTAL: per thread, the spikes sent and their digest (check 5)
  std::vector< std::uint64_t > experimental_spike_count_;
  std::vector< std::uint64_t > experimental_spike_digest_;

  // EXPERIMENTAL, NOT STOCK: the asynchronous exchange
  bool async_ = false;
  //! negative control (NEST_ASYNC_BREAK=reverse_ranks): the ranks' spikes delivered in reverse order
  bool async_break_reverse_ = false;
  //! checks only: the largest pseudo-random delay before a post, in microseconds (NEST_ASYNC_JITTER_US)
  long async_jitter_us_ = 0;
  //! checks only: measure how far other ranks ran ahead (NEST_ASYNC_STATS)
  bool async_stats_ = false;
  //! counts for the check: stamps posted, messages sent and received, the largest lead of another rank in slices
  size_t async_stamps_ = 0;
  size_t async_messages_sent_ = 0;
  size_t async_messages_received_ = 0;
  long async_largest_lead_ = -1;
  //! EXPERIMENTAL: slices that began with another rank's messages in for 0, 1 or 2 slices later (check 7)
  size_t async_lead_histogram_[ 3 ] = { 0, 0, 0 };
  //! the lockstep arm: the same messages, completed inside the slice
  bool lockstep_p2p_ = false;
  //! the master's cumulative exchange times in seconds, and the bytes it posted
  double async_time_post_ = 0.0;
  double async_time_wait_ = 0.0;
  double async_time_receive_ = 0.0;
  double async_time_complete_ = 0.0;
  size_t async_bytes_sent_ = 0;
  //! stamps whose messages are in flight
  std::set< long > async_expected_;
  //! stamps received by async_settle(), held for the slice that delivers them, with their chunk size
  std::map< long, std::pair< std::vector< SpikeData >, size_t > > async_held_;
#ifdef HAVE_MPI
  //! send buffers and requests of the messages posted, by stamp
  std::map< long, std::pair< std::vector< std::vector< SpikeData > >, std::vector< MPI_Request > > > async_sent_;
#endif
  //! receive the stamp's messages from every rank and build the buffer stock delivery reads
  void async_receive_( const long stamp );
  //! complete every send posted for stamps up to `stamp`
  void async_complete_sends_( const long stamp );
  std::vector< OffGridSpikeData > send_buffer_off_grid_spike_data_;
  std::vector< OffGridSpikeData > recv_buffer_off_grid_spike_data_;

  std::vector< TargetData > send_buffer_target_data_;
  std::vector< TargetData > recv_buffer_target_data_;

  //! whether size of MPI buffer for communication of connections was changed
  bool buffer_size_target_data_has_changed_;

  /**
   * Largest number of spikes sent from any rank to any other rank in last spike exchange round.
   *
   * The spike buffer section for any rank must be at least this size. Therefore, this number controls
   * buffer resizing.
   */
  size_t global_max_spikes_per_rank_;

  double send_recv_buffer_shrink_limit_;  //!< shrink buffer only if below this limit
  double send_recv_buffer_shrink_spare_;  //!< leave this fraction more space than minimally needed
  double send_recv_buffer_grow_extra_;    //!< when growing, add this fraction extra space

  /**
   * Log all resize events.
   *
   * This is maintained by the main thread, which is responsible for communication and resizing.
   */
  BufferResizeLog send_recv_buffer_resize_log_;

  PerThreadBoolIndicator gather_completed_checker_;

  // private stop watches for benchmarking purposes
  // (intended for internal core developers, not for use in the public API)
  Stopwatch< StopwatchGranularity::Detailed, StopwatchParallelism::MasterOnly > sw_collocate_spike_data_;
  Stopwatch< StopwatchGranularity::Detailed, StopwatchParallelism::MasterOnly > sw_communicate_spike_data_;
  Stopwatch< StopwatchGranularity::Detailed, StopwatchParallelism::MasterOnly > sw_communicate_target_data_;
};

inline void
EventDeliveryManager::reset_spike_register_( const size_t tid )
{
  emitted_spikes_register_[ tid ]->clear();
  off_grid_emitted_spikes_register_[ tid ]->clear();
}

inline bool
EventDeliveryManager::is_marked_for_removal_( const Target& target )
{
  return target.is_processed();
}

inline void
EventDeliveryManager::send_to_node( Event& e )
{
  e();
}

inline bool
EventDeliveryManager::get_off_grid_communication() const
{
  return off_grid_spiking_;
}

inline void
EventDeliveryManager::set_off_grid_communication( bool off_grid_spiking )
{
  off_grid_spiking_ = off_grid_spiking;
}

inline size_t
EventDeliveryManager::read_toggle() const
{
  // define in terms of write_toggle() to ensure consistency
  return 1 - write_toggle();
}

inline long
EventDeliveryManager::get_modulo( long d )
{
  // Note, here d may be 0, since bin 0 represents the "current" time
  // when all events due are read out.
  assert( static_cast< std::vector< long >::size_type >( d ) < moduli_.size() );

  return moduli_[ d ];
}

inline long
EventDeliveryManager::get_slice_modulo( long d )
{
  // Note, here d may be 0, since bin 0 represents the "current" time
  // when all events due are read out.
  assert( static_cast< std::vector< long >::size_type >( d ) < slice_moduli_.size() );

  return slice_moduli_[ d ];
}

}  // namespace nest

#endif /* EVENT_DELIVERY_MANAGER_H */
