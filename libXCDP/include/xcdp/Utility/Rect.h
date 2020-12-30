#pragma once

#include <vector>

#include "xcdp/Utility/Interval.h"
#include "xcdp/Utility/vec2.h"

namespace xcdp {

struct Rect 
	{
	Rect();
	Rect( const Interval & d, const Interval & r );
	Rect( float x1, float y1, float x2, float y2 );
	Rect( const std::vector<vec2> & pointsToFit, float expansion = 1.2 );

	bool operator!=( const Rect & o ) const;
	bool valid() const;
	Rect intersect( const Rect & other ) const;
	bool contains( const vec2 & p ) const;
	bool contains( float x, float y ) const;
	float aspect() const;

	float x1() const;
	float x2() const;
	float y1() const;
	float y2() const;
	float w() const;
	float h() const;

	Interval d, r;

	static const Rect Empty;
	static const Rect R2;
	};

inline const Rect Rect::Empty( Interval::Empty, Interval::Empty );
inline const Rect Rect::R2( Interval::R, Interval::R );

}