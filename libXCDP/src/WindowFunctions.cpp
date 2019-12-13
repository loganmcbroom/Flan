#include "WindowFunctions.h"
#include <cmath>

#include "RealFunc.h"

const double pi = std::acos( -1.0 );

namespace xcdp::window {

double Hann( double x )
	{
	return 0.5 * ( 1.0 - cos( 2.0 * pi * x ) );
	}

}