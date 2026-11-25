/*
 *  order_certified_input.h
 *
 *  EXPERIMENTAL, NOT STOCK NEST 3.10 (branch EXPERIMENTAL-NOT-STOCK-order-certified-step).
 *
 *  Restructuring U2 of the order-certificate project: during a run, the input slots of a thread's iaf_psc_exp neurons
 *  live in one time-major table, [slot][neuron in node order], each entry the same std::array< double, 4 > as the
 *  neuron's own MultiChannelInputBuffer slot. The update reads slot t of all its neurons as one stream.
 *
 *  Why the results stay bit-identical: this is a change of layout only. Every slot receives the same additions in the
 *  same order, through the same statements, and is read and reset by the same statements (THEORY.md: layout moves keep
 *  every order).
 *
 *  Lifetime: the table is built from the neurons' own buffers at the first slice of each run, and written back to them,
 *  then released, at the end of the run. Between runs all input lives in stock's structures, so deliveries at the
 *  start of a run, a kernel reset or new nodes never meet a stale table.
 */

#ifndef ORDER_CERTIFIED_INPUT_H
#define ORDER_CERTIFIED_INPUT_H

// C++ includes:
#include <array>
#include <cstddef>
#include <vector>

namespace nest
{
class iaf_psc_exp;

namespace order_certified_input
{

struct ThreadInput
{
  bool tried = false;  //!< a build was attempted in this run (it is not retried every slice)
  bool built = false;
  size_t n = 0;                                //!< neurons in the table
  size_t n_slots = 0;                          //!< slots per neuron, as in the neurons' own buffers
  std::vector< iaf_psc_exp* > neurons;         //!< in the thread's node order
  std::vector< std::array< double, 4 > > data;  //!< [slot * n + neuron]
};

//! Size the tables once and install the end-of-run hook. Called from pre_run_hook.
void register_hooks();

//! The table of thread tid.
ThreadInput& table( size_t tid );

//! Build thread tid's table if it is not built (called by the update-phase prepass).
void prepass( size_t tid );

//! One slot of one neuron in the table.
inline std::array< double, 4 >&
slot( ThreadInput& t, size_t s, size_t i )
{
  return t.data[ s * t.n + i ];
}

}  // namespace order_certified_input
}  // namespace nest

#endif  // ORDER_CERTIFIED_INPUT_H
