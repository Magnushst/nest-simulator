/*
 *  cycle_time_log.h
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

#ifndef CYCLE_TIME_LOG_H
#define CYCLE_TIME_LOG_H

// C++ includes:
#include <vector>

// Includes from libnestutil:
#include "dictionary.h"

namespace nest
{

/**
 * Collect cycle duration and spike-count information for individual update cycles.
 */
class CycleTimeLog
{
public:
  CycleTimeLog();
  void clear();
  void add_entry( double cycle_update_time,
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
    long spike_buffer_last );
  //! The OpenMP synchronisation values of the same cycle, added once the threads' stopped values are published.
  void add_omp_entry( double omp_max, double omp_mean );
  //! EXPERIMENTAL, NOT STOCK: the master's point-to-point exchange in the same cycle (zero in NEST's own path)
  void add_exchange_entry( double post, double wait, double receive, double complete, long messages, long bytes );
  //! EXPERIMENTAL: whether cycle `cycle` lies in the per-thread window (NEST_CYCLE_THREAD_WINDOW=start:count)
  bool in_thread_window( size_t cycle ) const;
  //! EXPERIMENTAL: append each thread's time in one phase (0 update, 1 delivery, 2 OpenMP synchronisation)
  void add_thread_values( int phase, const std::vector< double >& end, const std::vector< double >& start );
  //! the number of cycles logged so far
  size_t
  num_cycles() const
  {
    return cycle_update_time_.size();
  }
  void to_dict( Dictionary& ) const;

private:
  std::vector< double > cycle_update_time_;    //!< Time of one update cycle
  std::vector< double > communicate_time_;     // Time of communicate in current cycle
  std::vector< size_t > local_spike_counter_;  // Local spike count in current cycle
  std::vector< double > update_max_;            // Per-cycle team maximum of neuron update
  std::vector< double > update_mean_;           // Per-cycle team mean of neuron update
  std::vector< double > deliver_max_;           // Per-cycle team maximum of spike delivery
  std::vector< double > deliver_mean_;          // Per-cycle team mean of spike delivery
  std::vector< double > omp_max_;               // Per-cycle team maximum of omp synchronisation
  std::vector< double > omp_mean_;              // Per-cycle team mean of omp synchronisation
  std::vector< double > gather_time_;           // Time of gather in current cycle
  std::vector< double > collocate_time_;        // Time of collocate in current cycle
  std::vector< double > mpi_time_;              // Time of mpi synchronisation in current cycle
  std::vector< size_t > deliver_connections_;   // Connections visited by spike delivery in current cycle
  std::vector< size_t > spike_buffer_first_;    // Spike buffer per rank in the first all-to-all of current cycle
  std::vector< size_t > spike_buffer_last_;     // Spike buffer per rank in the last all-to-all of current cycle
  // EXPERIMENTAL, NOT STOCK: the point-to-point exchange per cycle, and the per-thread window
  std::vector< double > exchange_post_;
  std::vector< double > exchange_wait_;
  std::vector< double > exchange_receive_;
  std::vector< double > exchange_complete_;
  std::vector< size_t > exchange_messages_;
  std::vector< size_t > exchange_bytes_;
  long thread_window_start_ = -1;
  long thread_window_count_ = 0;
  size_t thread_window_threads_ = 0;
  std::vector< double > thread_update_;
  std::vector< double > thread_deliver_;
  std::vector< double > thread_omp_;
};

}

#endif /* cycle_time_log_h */
