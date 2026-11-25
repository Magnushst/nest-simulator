/*
 *  order_certified_drive.h
 *
 *  EXPERIMENTAL, NOT STOCK NEST 3.10 (branch EXPERIMENTAL-NOT-STOCK-order-certified-step).
 *
 *  Restructuring U1 of the order-certificate project: the Poisson drive of the update phase is drawn at the start of
 *  the update phase, in the reference order (generators in the thread's node order, then lag, then connection order),
 *  and each target adds weight x count to its own input slot during its own update, instead of receiving one
 *  DSSpikeEvent per target from the generator's update.
 *
 *  Why the results stay bit-identical (the project's THEORY.md, Theorem 1):
 *  - within the update phase only devices write input slots, always for later slices, and no node reads them, so
 *    every slot receives the same additions in the same order;
 *  - the thread's random stream gives the same draws to the same targets, because only the generators draw in the
 *    update phase on an eligible thread (build() checks it);
 *  - the product weight x count is kept apart from the addition, as in the stock binary (no fused multiply-add).
 *
 *  A thread that fails any condition of build() keeps the stock path unchanged.
 */

#ifndef ORDER_CERTIFIED_DRIVE_H
#define ORDER_CERTIFIED_DRIVE_H

// C++ includes:
#include <cstddef>
#include <vector>

// Includes from nestkernel:
#include "nest_time.h"

namespace nest
{
class poisson_generator;

namespace order_certified_drive
{

//! The drive of one thread: its generators in node order and, per target, the connection's weight and delay.
struct ThreadTable
{
  bool dirty = true;     //!< rebuilt at the next update phase (set by pre_run_hook of generators and neurons)
  bool enabled = false;  //!< true if this thread's drive is drawn by the prepass
  std::vector< poisson_generator* > generators;      //!< in the thread's node order
  std::vector< std::vector< size_t > > gen_targets;  //!< per generator: target indices, in connection order
  std::vector< double > weight;                      //!< per target index: the connection's weight
  std::vector< long > delay;                         //!< per target index: the connection's delay in steps
  std::vector< unsigned long > counts;               //!< [target * stride + lag - from]: the draw; allocated once per build
  long stride = 0;                                   //!< counts per target row: the longest slice (min_delay steps)
};

//! Size the tables once (call_once) and install the prepass hook. Called from pre_run_hook.
void register_hook();

//! Mark the thread's table for rebuilding. Called from pre_run_hook.
void mark_dirty( size_t tid );

//! The table of thread tid; only valid after register_hook().
ThreadTable& table( size_t tid );

//! True if the prepass draws the drive of thread tid.
bool handles_thread( size_t tid );

//! The prepass: rebuild the table if needed, then draw this slice's counts in the reference order.
void prepass( size_t tid, const Time& origin, long from, long to );

/**
 * Keep a product out of any fused multiply-add: the value is materialised in a register, so the compiler cannot
 * contract it into the following addition. The stock binary computes weight x multiplicity unfused (check 2).
 */
inline double
materialise( double x )
{
#if defined( __x86_64__ )
  asm volatile( "" : "+x"( x ) );
#elif defined( __aarch64__ )
  asm volatile( "" : "+w"( x ) );
#else
  volatile double y = x;
  x = y;
#endif
  return x;
}

}  // namespace order_certified_drive
}  // namespace nest

#endif  // ORDER_CERTIFIED_DRIVE_H
