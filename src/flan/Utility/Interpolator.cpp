#include "flan/Utility/Interpolator.h"

#include <cmath>
#include <algorithm>

#include "spline/spline.h"

using namespace flan;

static const float pi = std::acos( -1.0f );
static const float sqrt2 = std::sqrt( 2.0f );

/** Constantly 0.5. */
const Interpolator Interpolators::midpoint = []( float x ) 
	{ 
	return 0.5f;
	};

/** Returns nearest integer. */
const Interpolator Interpolators::nearest = []( float x ) 
	{ 
	return std::round( x );
	};

/** Constantly 0.0. */
const Interpolator Interpolators::floor = []( float x ) 
	{ 
	return 0.0f;
	};

/** Constantly 1.0. */
const Interpolator Interpolators::ceil = []( float x ) 
	{ 
	return 1.0f;
	};

/** Input returning function. */
const Interpolator Interpolators::linear = []( float x ) 
	{ 
	return x;
	};

/** Smoothstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
const Interpolator Interpolators::smoothstep = []( float x ) 
	{ 
	return x * x * ( 3.0f - 2.0f * x );
	};

/** Smootherstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
const Interpolator Interpolators::smootherstep = []( float x ) 
	{ 
	return x * x * x * ( x * ( x * 6.0f - 15.0f ) + 10.0f );
	};

/** Sine based interpolation. This moves from a minimum to the following maximum of a sinusoid. */
const Interpolator Interpolators::sine = []( float x ) 
	{ 
	return ( 1.0f - std::cos( pi * x ) ) / 2.0f; // Don't tell anyone, I used cosine for the sine interpolator >=)
	};

const Interpolator Interpolators::sine2 = []( float x )
	{
	return sqrt2 * sin( pi / 4.0f * x );
	};

/** Square root */
const Interpolator Interpolators::sqrt = []( float x ) 
	{ 
	return std::sqrt( x ); 
	};

Function<float, float> flan::interpolate_points( const std::vector< std::pair< float, float > > & ps, Interpolator interp )
	{
	return [ps, interp]( float t )
		{
		if( ps.size() == 0 ) return 0.0f;
		if( t <= ps.front().first ) return ps.front().second;
		if( ps.back().first <= t ) return ps.back().second;

		const auto p2 = std::lower_bound( ps.begin(), ps.end(), t, []( const std::pair< float, float > & p, float t ){ return p.first < t; } );
		const auto p1 = p2 - 1;
		const float mix = interp( ( t - p1->first ) / ( p2->first - p1->first ) );
		return ( 1.0f - mix ) * p1->second + mix * p2->second;
		};
	}

Function<float, float> flan::interpolate_intervals( float deltaX, const std::vector< float > ys, Interpolator interp )
	{
	std::vector<std::pair< float, float >> points( ys.size() );

	float xAccum = 0;
	std::transform( ys.begin(), ys.end(), points.begin(), [deltaX, &xAccum]( float y )
		{ 
		auto point = std::pair<float,float>( xAccum, y );
		xAccum += deltaX;
		return point;
		} );

	return interpolate_points( points, interp );
	}

Function<float, float> flan::spline( const std::vector< std::pair< float, float > > points )
	{
	// Split points into xs and ys
	std::vector<double> xs( points.size() );
	std::vector<double> ys( points.size() );
	std::transform( points.begin(), points.end(), xs.begin(), []( std::pair<float,float> p ){ return p.first; } );
	std::transform( points.begin(), points.end(), ys.begin(), []( std::pair<float,float> p ){ return p.second; } );

	// Generate tk::spline from xs and ys
	tk::spline spline;
	spline.set_points( xs, ys );

	return [spline = std::move( spline )]( float t )
		{
		return spline( t );
		};
	}