#pragma once

#include <functional>

namespace flan {

/** Interpolators describe how values between known data points should approximate those points near it.
 * These functions are passed values on [0,1] and are expected to return values on that same range.
 * OpenCL based algorithms sample Iterators passed to them, so rapidly changing user-defined
 * Iterators may not work as expected with those algorithms.
 */
using Interpolator = std::function< float ( float ) >;

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

	/** Sine based interpolation. This moves from a minimum to the following maximum of a sinusoid. */
	extern const Interpolator sine;

	/** Square root */
	extern const Interpolator sqrt;
	}

}