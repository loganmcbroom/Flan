#include "flan/Utility/Interpolator.h"

#include <cmath>
#include <algorithm>

#include "spline/spline.h"

using namespace flan;

static const float pi = std::acos( -1.0f );
static const float sqrt2 = std::sqrt( 2.0f );

/** Constantly 0.5. */
Interpolator Interpolator::midpoint()
    {
    return []( float x ) 
		{ 
		return 0.5f;
		};
	}

/** Returns nearest integer. */
Interpolator Interpolator::nearest()
    {
    return []( float x ) 
		{ 
		return std::round( x );
		};
	}

/** Constantly 0.0. */
Interpolator Interpolator::floor()
    {
    return []( float x ) 
		{ 
		return 0.0f;
		};
	}

/** Constantly 1.0. */
Interpolator Interpolator::ceil()
    {
    return []( float x ) 
		{ 
		return 1.0f;
		};
	}

/** Input returning function. */
Interpolator Interpolator::linear()
    {
    return []( float x ) 
		{ 
		return x;
		};
	}

/** Smoothstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
Interpolator Interpolator::smoothstep()
    {
    return []( float x ) 
		{ 
		return x * x * ( 3.0f - 2.0f * x );
		};
	}

/** Smootherstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
Interpolator Interpolator::smootherstep()
    {
    return []( float x ) 
		{ 
		return x * x * x * ( x * ( x * 6.0f - 15.0f ) + 10.0f );
		};
	}

/** Sine based interpolation. This moves from a minimum to the following maximum of a sinusoid. */
Interpolator Interpolator::sine()
    {
    return []( float x ) 
		{ 
		return ( 1.0f - std::cos( pi * x ) ) / 2.0f; // Don't tell anyone, I used cosine for the sine interpolator >=)
		};
	}

Interpolator Interpolator::sine2()
    {
    return []( float x )
		{
		return sqrt2 * sin( pi / 4.0f * x );
		};
	}

/** Square root */
Interpolator Interpolator::sqrt()
    {
    return []( float x ) 
		{ 
		return std::sqrt( x ); 
		};
	}

Function<float, float> flan::interpolate_points( const std::vector<vec2> & ps, Interpolator && interp )
	{
	auto policy = interp.f.get_execution_policy();
	// What's going on, eh? Interp in non-copyable, so trying to move it into the lambda capture would create a non-copyable lambda.
	// Function uses std::function which must call the copy ctor of the callable object it holds, which it can't do! So we need to 
	// grant interp copyability.
	auto p_interp = std::make_shared<Interpolator>( std::move( interp ) );
	return Function<float, float>( [ps, p_interp]( float t )
		{
		if( ps.size() == 0 ) return 0.0f;
		if( t <= ps.front().x() ) return ps.front().y();
		if( ps.back().x() <= t ) return ps.back().y();

		const auto p2 = std::lower_bound( ps.begin(), ps.end(), t, []( const vec2 & p, float t ){ return p.x() < t; } );
		const auto p1 = p2 - 1;
		const float mix = (*p_interp)( ( t - p1->x() ) / ( p2->x() - p1->x() ) );
		return ( 1.0f - mix ) * p1->y() + mix * p2->y();
		}, policy );
	}

Function<float, float> flan::interpolate_intervals( float deltaX, const std::vector< float > ys, Interpolator && interp )
	{
	std::vector<vec2> points( ys.size() );

	float xAccum = 0;
	std::transform( ys.begin(), ys.end(), points.begin(), [deltaX, &xAccum]( float y )
		{ 
		auto point = vec2( xAccum, y );
		xAccum += deltaX;
		return point;
		} );

	return interpolate_points( points, std::move( interp ) );
	}

Function<float, float> flan::spline( const std::vector<vec2> points )
	{
	// Split points into xs and ys
	std::vector<double> xs( points.size() );
	std::vector<double> ys( points.size() );
	std::transform( points.begin(), points.end(), xs.begin(), []( vec2 p ){ return p.x(); } );
	std::transform( points.begin(), points.end(), ys.begin(), []( vec2 p ){ return p.y(); } );

	// Generate tk::spline from xs and ys
	tk::spline spline;
	spline.set_points( xs, ys );

	return [spline = std::move( spline )]( float t )
		{
		return spline( t );
		};
	}