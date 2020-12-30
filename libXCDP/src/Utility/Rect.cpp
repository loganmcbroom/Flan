#include "xcdp/Utility/Rect.h"

#include <algorithm>

namespace xcdp {

Rect::Rect() 
	: d( 0, 1 )
	, r( 0, 1 )
	{}

Rect::Rect( const Interval & _d, const Interval & _r )
	: d( _d )
	, r( _r )
	{}

Rect::Rect( float x1, float y1, float x2, float y2 )
	: d( x1, x2 )
	, r( y1, y2 )
	{}

Rect::Rect( const std::vector<vec2> & ps, float expansion )
	: d( 0, 0 )
	, r( 0, 0 )
	{
	auto xComp = []( const vec2 & p,  const vec2 & s ){ return p.x() < s.x(); }; 
	auto yComp = []( const vec2 & p,  const vec2 & s ){ return p.y() < s.y(); }; 

	const float left	= std::min_element( ps.begin(), ps.end(), xComp )->x();
	const float right	= std::max_element( ps.begin(), ps.end(), xComp )->x();
	const float bottom	= std::min_element( ps.begin(), ps.end(), yComp )->y();
	const float top		= std::max_element( ps.begin(), ps.end(), yComp )->y();

	const float xCenter = ( left + right ) / 2.0f;
	const float yCenter = ( top + bottom ) / 2.0f;

	d = {
		left	* expansion + xCenter * ( 1.0f - expansion ),
		right	* expansion + xCenter * ( 1.0f - expansion )
		};
	r = {
		bottom	* expansion + yCenter * ( 1.0f - expansion ),
		top		* expansion + yCenter * ( 1.0f - expansion ) 
		};
	
	}

bool Rect::operator!=( const Rect & o ) const { return o.d != d || o.r != r; }

bool Rect::valid() const { return d.valid() && r.valid(); }

Rect Rect::intersect( const Rect & other ) const
	{
	return { d.intersect( other.d ), r.intersect( other.r ) };
	}

bool Rect::contains( const vec2 & p ) const { return d.contains( p.x() ) && r.contains( p.y() ); }
bool Rect::contains( float x, float y ) const { return contains({x,y}); }

float Rect::aspect() const { return w() / h(); }

float Rect::x1() const { return d.x1; }
float Rect::x2() const { return d.x2; }
float Rect::y1() const { return r.x1; }
float Rect::y2() const { return r.x2; }
float Rect::w() const { return d.w(); }
float Rect::h() const { return r.w(); }

}