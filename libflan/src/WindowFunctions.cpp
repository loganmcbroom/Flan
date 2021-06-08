#include "flan/WindowFunctions.h"
#include <cmath>

#include "flan/Function.h"

namespace flan::Windows {

const float pi = std::acos( -1.0f );

float Hann( float x )
	{
	return 0.5f * ( 1.0f - cos( 2.0f * pi * x ) );
	}

float HannDFT2( float f )
	{
	if( f == 0 ) return 1.0f;
	if( std::abs( f ) == 1.0f ) return 0.5f;
	return std::sin( pi * f ) / ( pi * f * ( 1.0f - f * f ) );
	};

}