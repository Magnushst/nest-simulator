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
 * counter (draw, then three words naming the connection), so a connection regenerates bit for bit whenever it is
 * delivered. Nothing here depends on the rest of NEST, so it is tested on its own.
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

// Includes from thirdparty:
#include "Random123/philox.h"

namespace nest
{

/**
 * The 32-bit words of one keyed stream: Philox4x64 blocks at counters (draw, c1, c2, c3), with the draw word counting
 * up from the counter the stream starts at.
 */
class ProceduralWords
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

  //! A uniform double on (0, 1] from two words: 53 bits plus one unit, so a logarithm never sees 0.
  double
  open_closed()
  {
    const std::uint64_t hi = next();
    return static_cast< double >( ( ( hi << 21 ) | ( next() >> 11 ) ) + 1 ) * 0x1p-53;
  }

  //! A uniform double on [0, 1) from two words.
  double
  closed_open()
  {
    const std::uint64_t hi = next();
    return static_cast< double >( ( hi << 21 ) | ( next() >> 11 ) ) * 0x1p-53;
  }

  //! A uniform integer on [0, n), n > 0, exactly: Lemire's multiply-shift on 64 bits, with rejection.
  std::uint64_t
  below( const std::uint64_t n )
  {
    __extension__ using u128 = unsigned __int128;
    const std::uint64_t threshold = ( 0 - n ) % n;
    for ( ;; )
    {
      const std::uint64_t hi = next();
      const u128 m = static_cast< u128 >( ( hi << 32 ) | next() ) * n;
      if ( static_cast< std::uint64_t >( m ) >= threshold )
      {
        return static_cast< std::uint64_t >( m >> 64 );
      }
    }
  }

private:
  Philox::key_type key_;
  Philox::ctr_type counter_;
  Philox::ctr_type block_ {};
  unsigned used_ = 8;
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

  double
  normal( ProceduralWords& w ) const
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

  double
  draw( const ProceduralZiggurat& z, ProceduralWords& w ) const
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
