#include "flan/Utility/Interval.h"

#include "flan/Utility/Rect.h"

namespace flan {

const Interval Interval::Empty( 0, 0 );
const Interval Interval::R( -std::numeric_limits<float>::max(), std::numeric_limits<float>::max() );

Interval::Interval( float _x1, float _x2 )
	: x1( _x1 )
	, x2( _x2 )
	{}

bool Interval::operator!=( const Interval & o ) const 
	{ 
	return o.x1 != x1 || o.x2 != x2; 
	}

bool Interval::valid() const 
	{ 
	return x1 <= x2; 
	}

float Interval::w() const 
	{ 
	return x2 - x1; 
	}

float Interval::midpoint() const 
	{ 
	return ( x1 + x2 ) / 2.0f; 
	}

Rect Interval::operator*( const Interval & o ) const 
	{ 
	return { x1, o.x1, x2, o.x2 }; 
	}

Interval Interval::intersect( const Interval & other ) const 
	{
	const Interval newI = { std::max( x1, other.x1 ), std::min( x2, other.x2 ) };
	return newI.valid() ? newI : Empty;
	}

bool Interval::contains( float x ) const 
	{ 
	return x1 <= x && x < x2; 
	}

}