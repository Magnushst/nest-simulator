/*
 *  procedural_random.h
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
 * Keyed random draws for regenerated connections. Every draw is a pure function of a key (seed, projection) and a
 * counter naming what it is for, so a connection regenerates bit for bit whenever it is delivered. A projection's
 * shares over its sources use Philox4x64 (ProceduralWords); each connection's target, weight and delay use
 * Philox4x32-10 (ProceduralWords32), whose 32 x 32-bit products let eight connections be drawn at once with AVX2
 * (procedural_philox_lanes). Nothing here depends on the rest of NEST, so it is tested on its own.
 */

#ifndef PROCEDURAL_RANDOM_H
#define PROCEDURAL_RANDOM_H

// C++ includes:
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

#if defined( __AVX2__ )
#include <immintrin.h>
#endif

// Includes from thirdparty:
#include "Random123/philox.h"

namespace nest
{

//! Uniform numbers from a stream of 32-bit words: a word stream derives from this and gives next().
template < class Words >
class ProceduralUniforms
{
public:
  //! A uniform double on (0, 1] from two words: 53 bits plus one unit, so a logarithm never sees 0.
  double
  open_closed()
  {
    const std::uint64_t hi = words_().next();
    return static_cast< double >( ( ( hi << 21 ) | ( words_().next() >> 11 ) ) + 1 ) * 0x1p-53;
  }

  //! A uniform double on [0, 1) from two words.
  double
  closed_open()
  {
    const std::uint64_t hi = words_().next();
    return static_cast< double >( ( hi << 21 ) | ( words_().next() >> 11 ) ) * 0x1p-53;
  }

  //! A uniform integer on [0, n), n > 0, exactly: Lemire's multiply-shift on 64 bits, with rejection.
  std::uint64_t
  below( const std::uint64_t n )
  {
    __extension__ using u128 = unsigned __int128;
    const std::uint64_t threshold = ( 0 - n ) % n;
    for ( ;; )
    {
      const std::uint64_t hi = words_().next();
      const u128 m = static_cast< u128 >( ( hi << 32 ) | words_().next() ) * n;
      if ( static_cast< std::uint64_t >( m ) >= threshold )
      {
        return static_cast< std::uint64_t >( m >> 64 );
      }
    }
  }

private:
  Words&
  words_()
  {
    return static_cast< Words& >( *this );
  }
};

/**
 * The 32-bit words of one keyed stream: Philox4x64 blocks at counters (draw, c1, c2, c3), with the draw word counting
 * up from the counter the stream starts at.
 */
class ProceduralWords : public ProceduralUniforms< ProceduralWords >
{
public:
  using Philox = r123::Philox4x64;

  ProceduralWords( const Philox::key_type& key, const Philox::ctr_type& counter )
    : key_( key )
    , counter_( counter )
  {
  }

  std::uint32_t
  next()
  {
    if ( used_ == 8 )
    {
      block_ = Philox()( counter_, key_ );
      ++counter_[ 0 ];
      used_ = 0;
    }
    const std::uint64_t word = block_[ used_ / 2 ];
    const std::uint32_t r =
      ( used_ & 1u ) ? static_cast< std::uint32_t >( word >> 32 ) : static_cast< std::uint32_t >( word );
    ++used_;
    return r;
  }

private:
  Philox::key_type key_;
  Philox::ctr_type counter_;
  Philox::ctr_type block_ {};
  unsigned used_ = 8;
};

using ProceduralCtr32 = std::array< std::uint32_t, 4 >;
using ProceduralKey32 = std::array< std::uint32_t, 2 >;

/**
 * Philox4x32-10 (Salmon and others, SC 2011), as Random123's r123::Philox4x32 computes it (the tests compare the two)
 * and as the prototype engine writes it out, so that procedural_philox_lanes can give the same numbers eight at once.
 */
inline ProceduralCtr32
procedural_philox32( ProceduralCtr32 c, ProceduralKey32 k )
{
  constexpr std::uint64_t m0 = 0xD2511F53u;
  constexpr std::uint64_t m1 = 0xCD9E8D57u;
  for ( int r = 0; r < 10; ++r )
  {
    if ( r > 0 )
    {
      k[ 0 ] += 0x9E3779B9u;
      k[ 1 ] += 0xBB67AE85u;
    }
    const std::uint64_t p0 = m0 * c[ 0 ];
    const std::uint64_t p1 = m1 * c[ 2 ];
    c = { static_cast< std::uint32_t >( p1 >> 32 ) ^ c[ 1 ] ^ k[ 0 ],
      static_cast< std::uint32_t >( p1 ),
      static_cast< std::uint32_t >( p0 >> 32 ) ^ c[ 3 ] ^ k[ 1 ],
      static_cast< std::uint32_t >( p0 ) };
  }
  return c;
}

//! Eight Philox4x32-10 counters, one array per counter word, replaced in place by their blocks.
constexpr std::size_t procedural_lanes = 8;
struct ProceduralLanes
{
  alignas( 32 ) std::array< std::array< std::uint32_t, procedural_lanes >, 4 > w;
};

/**
 * procedural_philox32 for eight counters at once: with AVX2 one register per counter word, the even and odd lanes'
 * 32 x 32 -> 64-bit products by _mm256_mul_epu32, recombined by blends, as the prototype engine computes them;
 * otherwise a loop. The numbers are exactly procedural_philox32's.
 */
inline void
procedural_philox_lanes( ProceduralLanes& c, ProceduralKey32 k )
{
#if defined( __AVX2__ )
  const auto load = [ & ]( const std::size_t i )
  { return _mm256_load_si256( reinterpret_cast< const __m256i* >( c.w[ i ].data() ) ); };
  __m256i c0 = load( 0 );
  __m256i c1 = load( 1 );
  __m256i c2 = load( 2 );
  __m256i c3 = load( 3 );
  const __m256i m0 = _mm256_set1_epi32( static_cast< int >( 0xD2511F53u ) );
  const __m256i m1 = _mm256_set1_epi32( static_cast< int >( 0xCD9E8D57u ) );
  // the high and low 32 bits of m * x for each 32-bit lane of x
  const auto mul = []( const __m256i m, const __m256i x, __m256i& hi, __m256i& lo )
  {
    const __m256i even = _mm256_mul_epu32( m, x );
    const __m256i odd = _mm256_mul_epu32( m, _mm256_srli_epi64( x, 32 ) );
    hi = _mm256_blend_epi32( _mm256_srli_epi64( even, 32 ), odd, 0xAA );
    lo = _mm256_blend_epi32( even, _mm256_slli_epi64( odd, 32 ), 0xAA );
  };
  for ( int r = 0; r < 10; ++r )
  {
    if ( r > 0 )
    {
      k[ 0 ] += 0x9E3779B9u;
      k[ 1 ] += 0xBB67AE85u;
    }
    __m256i hi0, lo0, hi1, lo1;
    mul( m0, c0, hi0, lo0 );
    mul( m1, c2, hi1, lo1 );
    c0 = _mm256_xor_si256( _mm256_xor_si256( hi1, c1 ), _mm256_set1_epi32( static_cast< int >( k[ 0 ] ) ) );
    c1 = lo1;
    c2 = _mm256_xor_si256( _mm256_xor_si256( hi0, c3 ), _mm256_set1_epi32( static_cast< int >( k[ 1 ] ) ) );
    c3 = lo0;
  }
  const auto store = [ & ]( const std::size_t i, const __m256i v )
  { _mm256_store_si256( reinterpret_cast< __m256i* >( c.w[ i ].data() ), v ); };
  store( 0, c0 );
  store( 1, c1 );
  store( 2, c2 );
  store( 3, c3 );
#else
  for ( std::size_t l = 0; l < procedural_lanes; ++l )
  {
    const ProceduralCtr32 r = procedural_philox32( { c.w[ 0 ][ l ], c.w[ 1 ][ l ], c.w[ 2 ][ l ], c.w[ 3 ][ l ] }, k );
    for ( std::size_t i = 0; i < 4; ++i )
    {
      c.w[ i ][ l ] = r[ i ];
    }
  }
#endif
}

/**
 * The 32-bit words of one connection: the two Philox4x32-10 blocks drawn for it (blocks 0 and 1 of its counter, often
 * eight connections at once), then further blocks of the same counter, its last word counting on from 2.
 */
class ProceduralWords32 : public ProceduralUniforms< ProceduralWords32 >
{
public:
  ProceduralWords32( const ProceduralKey32& key, const ProceduralCtr32& counter, const std::array< std::uint32_t, 8 >& first )
    : key_( key )
    , counter_( counter )
    , buffer_( first )
  {
  }

  std::uint32_t
  next()
  {
    if ( used_ == available_ )
    {
      counter_[ 3 ] = block_++;
      const ProceduralCtr32 b = procedural_philox32( counter_, key_ );
      std::copy( b.begin(), b.end(), buffer_.begin() );
      available_ = 4;
      used_ = 0;
    }
    return buffer_[ used_++ ];
  }

private:
  ProceduralKey32 key_;
  ProceduralCtr32 counter_;
  std::array< std::uint32_t, 8 > buffer_;
  unsigned used_ = 0;
  unsigned available_ = 8;
  std::uint32_t block_ = 2;
};

/**
 * Normal deviates by the ziggurat method with 256 layers, in Doornik's form (ZIGNOR, 2005), as the prototype engine
 * draws them: the layer and the uniform come from separate bits. R is the base layer's edge and V each layer's area
 * (Marsaglia and Tsang, 2000); the recursion closes the top layer.
 */
class ProceduralZiggurat
{
public:
  ProceduralZiggurat()
  {
    const auto f = []( const double v ) { return std::exp( -0.5 * v * v ); };
    x_[ 0 ] = v_ / f( r_ );
    x_[ 1 ] = r_;
    for ( std::size_t i = 2; i < layers_; ++i )
    {
      x_[ i ] = std::sqrt( -2.0 * std::log( v_ / x_[ i - 1 ] + f( x_[ i - 1 ] ) ) );
    }
    x_[ layers_ ] = 0.0;
    for ( std::size_t i = 0; i < layers_; ++i )
    {
      ratio_[ i ] = x_[ i + 1 ] / x_[ i ];
    }
  }

  template < class Words >
  double
  normal( Words& w ) const
  {
    for ( ;; )
    {
      const std::uint32_t a = w.next();
      const std::uint32_t b = w.next();
      const std::uint32_t i = a & 0xFFu;
      const bool negative = ( a >> 8 ) & 1u;
      const double u = static_cast< double >( ( static_cast< std::uint64_t >( a >> 9 ) << 30 ) | ( b >> 2 ) ) * 0x1p-53;
      if ( u < ratio_[ i ] )
      {
        return negative ? -u * x_[ i ] : u * x_[ i ];
      }
      if ( i == 0 ) // the tail beyond R, by Marsaglia's method
      {
        double x = 0.0;
        double y = 0.0;
        do
        {
          x = std::log( w.open_closed() ) / r_;
          y = std::log( w.open_closed() );
        } while ( -2.0 * y < x * x );
        return negative ? x - r_ : r_ - x;
      }
      const double x = u * x_[ i ];
      const double f0 = std::exp( -0.5 * ( x_[ i ] * x_[ i ] - x * x ) );
      const double f1 = std::exp( -0.5 * ( x_[ i + 1 ] * x_[ i + 1 ] - x * x ) );
      if ( f1 + w.closed_open() * ( f0 - f1 ) < 1.0 )
      {
        return negative ? -x : x;
      }
    }
  }

private:
  static constexpr std::size_t layers_ = 256;
  static constexpr double r_ = 3.6541528853610088;
  static constexpr double v_ = 4.92867323399e-3;
  std::array< double, layers_ + 1 > x_ {};
  std::array< double, layers_ > ratio_ {};
};

/**
 * The counter of a keyed stream: the draw word first, then two words the stream defines, then the stream's tag in the
 * top byte of the last word above the virtual process, so that no two streams of one key share a counter.
 */
enum class ProceduralStream : std::uint64_t
{
  connection = 1, //!< one connection's target, weight and delay (Philox4x32-10, procedural_connection_counter)
  share = 2       //!< (tree node, 0): a split of a projection's connections over its sources
};

inline ProceduralWords::Philox::ctr_type
procedural_counter( const ProceduralStream stream, const std::uint64_t a, const std::uint64_t b, const std::uint64_t vp )
{
  return ProceduralWords::Philox::ctr_type { { 0, a, b, ( static_cast< std::uint64_t >( stream ) << 56 ) | vp } };
}

/**
 * The Philox4x32-10 counter of connection i of source s on virtual process vp: the source, the stream's tag in the top
 * byte above the virtual process (below 2^24), the connection, and the block. The caller keeps s below 2^32 and i
 * below 2^32.
 */
inline ProceduralCtr32
procedural_connection_counter( const std::uint32_t s, const std::uint32_t vp, const std::uint32_t i, const std::uint32_t block )
{
  return { s, ( static_cast< std::uint32_t >( ProceduralStream::connection ) << 24 ) | vp, i, block };
}

/**
 * log(k!) minus its Stirling approximation, for BTRS's acceptance test (Hoermann 1993): a table below 10, the series
 * above. As the prototype engine computes it.
 */
inline double
procedural_stirling_tail( const double k )
{
  static constexpr double tail[ 10 ] = { 0.08106146679532726,
    0.04134069595540929,
    0.02767792568499834,
    0.02079067210376509,
    0.01664469118982119,
    0.01387612882307075,
    0.01189670994589177,
    0.01041126526197209,
    0.009255462182712733,
    0.008330563433362871 };
  if ( k <= 9 )
  {
    return tail[ static_cast< int >( k ) ];
  }
  const double kp1sq = ( k + 1 ) * ( k + 1 );
  return ( 1.0 / 12 - ( 1.0 / 360 - 1.0 / 1260 / kp1sq ) / kp1sq ) / ( k + 1 );
}

/**
 * Binomial(n, p), exactly, as the prototype engine draws it: inversion while n min(p, 1 - p) is below 10, else
 * Hoermann's transformed rejection with squeeze (BTRS, 1993); p above 1/2 by symmetry.
 */
inline std::uint64_t
procedural_binomial( const std::uint64_t n, const double p, ProceduralWords& w )
{
  if ( n == 0 or p <= 0.0 )
  {
    return 0;
  }
  if ( p >= 1.0 )
  {
    return n;
  }
  if ( p > 0.5 )
  {
    return n - procedural_binomial( n, 1.0 - p, w );
  }
  const double nd = static_cast< double >( n );
  const double q = 1.0 - p;
  if ( nd * p < 10.0 )
  {
    const double r = p / q;
    const double p0 = std::exp( nd * std::log1p( -p ) );
    for ( ;; ) // x beyond the rounded sum of all probabilities (about 1e-15) is drawn again
    {
      const double x = w.closed_open();
      double pk = p0;
      double f = pk;
      std::uint64_t k = 0;
      while ( x >= f and k < n and pk > 0.0 )
      {
        pk *= r * static_cast< double >( n - k ) / static_cast< double >( k + 1 );
        ++k;
        f += pk;
      }
      if ( x < f )
      {
        return k;
      }
    }
  }
  const double spq = std::sqrt( nd * p * q );
  const double b = 1.15 + 2.53 * spq;
  const double a = -0.0873 + 0.0248 * b + 0.01 * p;
  const double c = nd * p + 0.5;
  const double vr = 0.92 - 4.2 / b;
  const double alpha = ( 2.83 + 5.1 / b ) * spq;
  const double m = std::floor( ( nd + 1 ) * p );
  for ( ;; )
  {
    const double u = w.closed_open() - 0.5;
    const double v = w.closed_open();
    const double us = 0.5 - std::abs( u );
    if ( us <= 0.0 )
    {
      continue;
    }
    const double k = std::floor( ( 2.0 * a / us + b ) * u + c );
    if ( k < 0.0 or k > nd )
    {
      continue;
    }
    if ( us >= 0.07 and v <= vr )
    {
      return static_cast< std::uint64_t >( k );
    }
    const double lv = std::log( v * alpha / ( a / ( us * us ) + b ) );
    const double bound = ( m + 0.5 ) * std::log( ( m + 1 ) / ( p / q * ( nd - m + 1 ) ) )
      + ( nd + 1 ) * std::log( ( nd - m + 1 ) / ( nd - k + 1 ) )
      + ( k + 0.5 ) * std::log( p / q * ( nd - k + 1 ) / ( k + 1 ) ) + procedural_stirling_tail( m )
      + procedural_stirling_tail( nd - m ) - procedural_stirling_tail( k ) - procedural_stirling_tail( nd - k );
    if ( lv <= bound )
    {
      return static_cast< std::uint64_t >( k );
    }
  }
}

/**
 * One weight or delay law of a regenerated projection, as NEST's parameter objects define it: a constant, a normal,
 * or a lognormal (given by its underlying normal, as std::lognormal_distribution), then left alone, raised to low
 * (NEST's math.max), lowered to high (math.min), or drawn again while below low or above high (math.redraw, which
 * gives up after max_redraws draws, as NEST does).
 */
struct ProceduralLaw
{
  enum class Base
  {
    constant,
    normal,
    lognormal
  };
  enum class Clip
  {
    none,
    max,
    min,
    redraw
  };

  Base base = Base::constant;
  Clip clip = Clip::none;
  double mean = 0.0; //!< the constant, or the (underlying) normal's mean
  double sigma = 0.0; //!< the (underlying) normal's standard deviation
  double low = -std::numeric_limits< double >::infinity();
  double high = std::numeric_limits< double >::infinity();
  std::size_t max_redraws = 1000;

  template < class Words >
  double
  draw( const ProceduralZiggurat& z, Words& w ) const
  {
    for ( std::size_t draws = 1;; ++draws )
    {
      double x = mean;
      if ( base != Base::constant )
      {
        x = mean + sigma * z.normal( w );
        if ( base == Base::lognormal )
        {
          x = std::exp( x );
        }
      }
      switch ( clip )
      {
      case Clip::none:
        return x;
      case Clip::max:
        return std::max( x, low );
      case Clip::min:
        return std::min( x, high );
      case Clip::redraw:
        if ( not( x < low or x > high ) )
        {
          return x;
        }
        if ( draws == max_redraws )
        {
          throw std::runtime_error( "Number of redraws exceeded limit" );
        }
        break;
      }
    }
  }
};

} // namespace nest

#endif /* #ifndef PROCEDURAL_RANDOM_H */
