#include "flan/Utility/Interpolator.h"

using namespace flan;

static const float pi = std::acos( -1.0f );

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
	return ( 1.0f - cos( pi * x ) ) / 2.0f; 
	};

/** Square root */
const Interpolator Interpolators::sqrt = []( float x ) 
	{ 
	return std::sqrt( x ); 
	};