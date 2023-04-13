#pragma once

#include "flan/defines.h"

namespace flan {

/** This is the data type stored at each buffer position. It represents a frequency, and the magnitude of that frequency.
*/
struct MF 
	{ 
	// MF( vec2 v ) : m( v.x() ), f( v.y() ) {} // Ctor from vec2
	Magnitude m;
	Frequency f; 
	};

}