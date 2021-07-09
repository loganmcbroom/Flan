#pragma once

namespace flan {

struct Rect;

/** Half-open interval **/
struct Interval
	{
	Interval( float x1, float x2 );

	bool operator!=( const Interval & o ) const;

	bool valid() const;

	float w() const;

	float midpoint() const;

	Rect operator*( const Interval & o ) const;

	Interval intersect( const Interval & other ) const;

	bool contains( float x ) const;

	float x1, x2;

	static const Interval Empty;
	static const Interval R;
	};

}