#pragma once

#include "flan/defines.h"

namespace flan {

using UnsignedPitch = float;

struct Pitch {
	Pitch() : p( 0 ), positive_frequency( true ) {}
	Pitch( UnsignedPitch _p, bool sign = true ) : p( _p ), positive_frequency( sign ) {}

	Pitch operator+( UnsignedPitch q ) { return Pitch( p + q, positive_frequency ); }
	Pitch operator-( UnsignedPitch q ) { return Pitch( p - q, positive_frequency ); }
	Pitch operator*( float scale ) { return Pitch( p * scale, positive_frequency ); }

	UnsignedPitch p; 

	// Negative frequencies are seemingly unavoidable, but you can't take a log of a negative, so we need to 
	// log2( abs( freq ) ) to get pitch while also tracking the sign of the input frequency. This solution is sub-par,
	// so it may be changed later. For one it messes up the alignment of MP. Ideally we could steal a bit of pitch to 
	// store this info, but I don't know of any way to do that.
	bool positive_frequency; 
};

/** This is the data type stored at each buffer position. It represents a frequency, and the magnitude of that frequency.
*/
struct MP 
	{ 
	// MF( vec2 v ) : m( v.x() ), f( v.y() ) {} // Ctor from vec2
	Magnitude m;
	Pitch p;
	};

}