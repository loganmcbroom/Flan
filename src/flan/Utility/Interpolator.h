#pragma once

#include "flan/Function.h"

namespace flan {

/** Interpolators describe how values between known data points should approximate those points near it.
 * These functions are passed values on [0,1] and are expected to return values on that same range.
 */
struct Interpolator {

	Interpolator( const Interpolator & ) 				= delete;
	Interpolator & operator=( const Interpolator & ) 	= delete;
	Interpolator( Interpolator && ) 					= default;
	Interpolator & operator=( Interpolator && ) 		= default;
	~Interpolator() 									= default;

	template<typename F>
	Interpolator( const F & _f ) : f( _f ) {}

	Interpolator( Function<float, float> && _f ) : f( std::move( _f ) ) {}
	//Interpolator( const Function<float, float> & _f ) : f( _f.copy() ) {}

	float operator()( float x ) const { return f(x); }

	/** Constantly 0.5. */
	static Interpolator midpoint();

	/** Returns nearest integer. */
	static Interpolator nearest();

	/** Constantly 0.0. */
	static Interpolator floor();

	/** Constantly 1.0. */
	static Interpolator ceil();

	/** Input returning function. */
	static Interpolator linear();

	/** Smoothstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
	static Interpolator smoothstep();

	/** Smootherstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
	static Interpolator smootherstep();

	/** Sine based interpolation. This moves from a minimum to the following maximum of a sinusoid, shifted and scaled to have a range of [0,1]. */
	static Interpolator sine();

	/** This starts at the midpoint of a sinusoid and moves halfway to a peak, shifted and scaled to have a range of [0,1]. */
	static Interpolator sine2();

	/** Square root */
	static Interpolator sqrt();

	Function<float, float> f;
};

/** Generate a function that passes through a given set of points.
*
* \param points Points that the generated function must pass through.
* \param interp How the generated function should move between points.
*/
Function<float, float> interpolate_points( const std::vector<vec2> & points, Interpolator && interp = Interpolator::linear() );

/** Generate a function that passes through a given set of points.
*
* \param deltaX Distance between x samples, starting at x = 0.
* \param points Y values that the generated function must pass through.
* \param interp How the generated function should move between points.
*/
Function<float, float> interpolate_intervals( float deltaX, const std::vector< float > ys, Interpolator && interp = Interpolator::linear() );

/** Generate a cubic spline that passes through a given set of points.
	*
	* \param points Points that the generated function must pass through.
	*/
Function<float, float> spline( const std::vector<vec2> points );

}