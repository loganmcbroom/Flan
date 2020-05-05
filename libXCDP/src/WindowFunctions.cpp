#include "WindowFunctions.h"
#include <cmath>

#include "RealFunc.h"

namespace xcdp::window {

const float pi = std::acos( -1.0f );

float Hann( float x )
	{
	return 0.5f * ( 1.0f - cos( 2.0f * pi * x ) );
	}

}