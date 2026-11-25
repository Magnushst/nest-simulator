/*
 *  cycle_time_log.cpp
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

#include "cycle_time_log.h"

#include <utility>

#include <string>

#include <cstdlib>

#include <algorithm>

// Includes from nestkernel:
#include "kernel_manager.h"
#include "nest_names.h"

namespace nest
{

CycleTimeLog::CycleTimeLog()
  : cycle_update_time_()
  , communicate_time_()
  , local_spike_counter_()
  , update_max_()
  , update_mean_()
  , deliver_max_()
  , deliver_mean_()
  , omp_max_()
  , omp_mean_()
  , gather_time_()
  , collocate_time_()
  , mpi_time_()
  , deliver_connections_()
  , spike_buffer_first_()
  , spike_buffer_last_()
{
}

void
CycleTimeLog::clear()
{
  cycle_update_time_.clear();
  communicate_time_.clear();
  local_spike_counter_.clear();
  update_max_.clear();
  update_mean_.clear();
  deliver_max_.clear();
  deliver_mean_.clear();
  omp_max_.clear();
  omp_mean_.clear();
  gather_time_.clear();
  collocate_time_.clear();
  mpi_time_.clear();
  deliver_connections_.clear();
  spike_buffer_first_.clear();
  spike_buffer_last_.clear();
  // EXPERIMENTAL, NOT STOCK
  exchange_post_.clear();
  exchange_wait_.clear();
  exchange_receive_.clear();
  exchange_complete_.clear();
  exchange_messages_.clear();
  exchange_bytes_.clear();
  thread_update_.clear();
  thread_deliver_.clear();
  thread_omp_.clear();
  thread_window_start_ = -1;
  thread_window_count_ = 0;
  thread_window_threads_ = 0;
  const char* window = std::getenv( "NEST_CYCLE_THREAD_WINDOW" );
  if ( window != nullptr )
  {
    const std::string spec( window );
    const size_t colon = spec.find( ':' );
    if ( colon != std::string::npos )
    {
      thread_window_start_ = std::max( 0L, std::atol( spec.substr( 0, colon ).c_str() ) );
      thread_window_count_ = std::max( 0L, std::atol( spec.substr( colon + 1 ).c_str() ) );
    }
  }
}

void
CycleTimeLog::add_entry( double cycle_update_time,
  double communicate_time,
  long local_spike_counter,
  double update_max,
  double update_mean,
  double deliver_max,
  double deliver_mean,
  double gather_time,
  double collocate_time,
  double mpi_time,
  long deliver_connections,
  long spike_buffer_first,
  long spike_buffer_last )
{
  cycle_update_time_.emplace_back( cycle_update_time );
  communicate_time_.emplace_back( communicate_time );
  local_spike_counter_.emplace_back( local_spike_counter );
  update_max_.emplace_back( update_max );
  update_mean_.emplace_back( update_mean );
  deliver_max_.emplace_back( deliver_max );
  deliver_mean_.emplace_back( deliver_mean );
  gather_time_.emplace_back( gather_time );
  collocate_time_.emplace_back( collocate_time );
  mpi_time_.emplace_back( mpi_time );
  deliver_connections_.emplace_back( deliver_connections );
  spike_buffer_first_.emplace_back( spike_buffer_first );
  spike_buffer_last_.emplace_back( spike_buffer_last );
}

void
CycleTimeLog::add_exchange_entry( double post, double wait, double receive, double complete, long messages,
  long bytes )
{
  exchange_post_.emplace_back( post );
  exchange_wait_.emplace_back( wait );
  exchange_receive_.emplace_back( receive );
  exchange_complete_.emplace_back( complete );
  exchange_messages_.emplace_back( messages );
  exchange_bytes_.emplace_back( bytes );
}

bool
CycleTimeLog::in_thread_window( size_t cycle ) const
{
  return thread_window_start_ >= 0 and static_cast< long >( cycle ) >= thread_window_start_
    and static_cast< long >( cycle ) < thread_window_start_ + thread_window_count_;
}

void
CycleTimeLog::add_thread_values( int phase, const std::vector< double >& end, const std::vector< double >& start )
{
  thread_window_threads_ = end.size();
  std::vector< double >& out = phase == 0 ? thread_update_ : ( phase == 1 ? thread_deliver_ : thread_omp_ );
  for ( size_t t = 0; t < end.size(); ++t )
  {
    out.emplace_back( end[ t ] - start[ t ] );
  }
}

void
CycleTimeLog::add_omp_entry( double omp_max, double omp_mean )
{
  omp_max_.emplace_back( omp_max );
  omp_mean_.emplace_back( omp_mean );
}

void
CycleTimeLog::to_dict( Dictionary& events ) const
{
  auto& times = events.get_or_create_vector< double >( names::time_update );
  times.insert( times.end(), cycle_update_time_.begin(), cycle_update_time_.end() );

  auto& communicate_time = events.get_or_create_vector< double >( names::time_communicate_spike_data );
  communicate_time.insert( communicate_time.end(), communicate_time_.begin(), communicate_time_.end() );

  auto& local_spike_counter = events.get_or_create_vector< long >( names::local_spike_counter );
  local_spike_counter.insert( local_spike_counter.end(), local_spike_counter_.begin(), local_spike_counter_.end() );

  auto& update_max = events.get_or_create_vector< double >( names::phaselok_update_max );
  update_max.insert( update_max.end(), update_max_.begin(), update_max_.end() );

  auto& update_mean = events.get_or_create_vector< double >( names::phaselok_update_mean );
  update_mean.insert( update_mean.end(), update_mean_.begin(), update_mean_.end() );

  auto& deliver_max = events.get_or_create_vector< double >( names::phaselok_deliver_max );
  deliver_max.insert( deliver_max.end(), deliver_max_.begin(), deliver_max_.end() );

  auto& deliver_mean = events.get_or_create_vector< double >( names::phaselok_deliver_mean );
  deliver_mean.insert( deliver_mean.end(), deliver_mean_.begin(), deliver_mean_.end() );

  auto& omp_max = events.get_or_create_vector< double >( names::phaselok_omp_max );
  omp_max.insert( omp_max.end(), omp_max_.begin(), omp_max_.end() );

  auto& omp_mean = events.get_or_create_vector< double >( names::phaselok_omp_mean );
  omp_mean.insert( omp_mean.end(), omp_mean_.begin(), omp_mean_.end() );

  auto& gather_time = events.get_or_create_vector< double >( names::phaselok_gather );
  gather_time.insert( gather_time.end(), gather_time_.begin(), gather_time_.end() );

  auto& collocate_time = events.get_or_create_vector< double >( names::phaselok_collocate );
  collocate_time.insert( collocate_time.end(), collocate_time_.begin(), collocate_time_.end() );

  auto& mpi_time = events.get_or_create_vector< double >( names::phaselok_mpi );
  mpi_time.insert( mpi_time.end(), mpi_time_.begin(), mpi_time_.end() );

  auto& deliver_connections = events.get_or_create_vector< long >( names::phaselok_deliver_connections );
  deliver_connections.insert(
    deliver_connections.end(), deliver_connections_.begin(), deliver_connections_.end() );

  auto& spike_buffer_first = events.get_or_create_vector< long >( names::phaselok_spike_buffer_first );
  spike_buffer_first.insert( spike_buffer_first.end(), spike_buffer_first_.begin(), spike_buffer_first_.end() );

  auto& spike_buffer_last = events.get_or_create_vector< long >( names::phaselok_spike_buffer_last );
  spike_buffer_last.insert( spike_buffer_last.end(), spike_buffer_last_.begin(), spike_buffer_last_.end() );

  // EXPERIMENTAL, NOT STOCK: the point-to-point exchange per cycle, and the per-thread window
  const std::vector< std::pair< const std::string*, const std::vector< double >* > > doubles = {
    { &names::phaselok_exchange_post, &exchange_post_ },
    { &names::phaselok_exchange_wait, &exchange_wait_ },
    { &names::phaselok_exchange_receive, &exchange_receive_ },
    { &names::phaselok_exchange_complete, &exchange_complete_ },
    { &names::phaselok_thread_update, &thread_update_ },
    { &names::phaselok_thread_deliver, &thread_deliver_ },
    { &names::phaselok_thread_omp, &thread_omp_ } };
  for ( const auto& [ name, values ] : doubles )
  {
    auto& v = events.get_or_create_vector< double >( *name );
    v.insert( v.end(), values->begin(), values->end() );
  }
  auto& messages = events.get_or_create_vector< long >( names::phaselok_exchange_messages );
  messages.insert( messages.end(), exchange_messages_.begin(), exchange_messages_.end() );
  auto& bytes = events.get_or_create_vector< long >( names::phaselok_exchange_bytes );
  bytes.insert( bytes.end(), exchange_bytes_.begin(), exchange_bytes_.end() );
  auto& window = events.get_or_create_vector< long >( names::phaselok_thread_window );
  window.insert(
    window.end(), { thread_window_start_, thread_window_count_, static_cast< long >( thread_window_threads_ ) } );
}

}
