#include "flan/WindowFunctions.h"
#include <cmath>

#include "flan/Function.h"

namespace flan::Windows {

const float pi = std::acos( -1.0f );

float hann( float x )
	{
	return 0.5f * ( 1.0f - cos( 2.0f * pi * x ) );
	}

}