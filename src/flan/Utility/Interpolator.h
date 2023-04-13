#pragma once

#include <functional>

#include "flan/Function.h"

namespace flan {

/** Interpolators describe how values between known data points should approximate those points near it.
 * These functions are passed values on [0,1] and are expected to return values on that same range.
 * OpenCL based algorithms sample Iterators passed to them, so rapidly changing user-defined
 * Iterators may not work as expected with those algorithms.
 */
using Interpolator = std::function<float (float)>;

namespace Interpolators 
	{
	/** Constantly 0.5. */
	extern const Interpolator midpoint;

	/** Returns nearest integer. */
	extern const Interpolator nearest;

	/** Constantly 0.0. */
	extern const Interpolator floor;

	/** Constantly 1.0. */
	extern const Interpolator ceil;

	/** Input returning function. */
	extern const Interpolator linear;

	/** Smoothstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
	extern const Interpolator smoothstep;

	/** Smootherstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
	extern const Interpolator smootherstep;

	/** Sine based interpolation. This moves from a minimum to the following maximum of a sinusoid, shifted and scaled to have a range of [0,1]. */
	extern const Interpolator sine;

	/** This starts at the midpoint of a sinusoid and moves halfway to a peak, shifted and scaled to have a range of [0,1]. */
	extern const Interpolator sine2;

	/** Square root */
	extern const Interpolator sqrt;
	}

/** Generate a function that passes through a given set of points.
*
* \param points Points that the generated function must pass through.
* \param interp How the generated function should move between points.
*/
Func1x1 interpolatePoints( const std::vector< std::pair< float, float > > & points, Interpolator interp = Interpolators::linear );

/** Generate a function that passes through a given set of points.
*
* \param deltaX Distance between x samples, starting at x = 0.
* \param points Y values that the generated function must pass through.
* \param interp How the generated function should move between points.
*/
Func1x1 interpolateIntervals( float deltaX, const std::vector< float > ys, Interpolator interp = Interpolators::linear );

/** Generate a cubic spline that passes through a given set of points.
	*
	* \param points Points that the generated function must pass through.
	*/
Func1x1 spline( const std::vector< std::pair< float, float > > points );

}