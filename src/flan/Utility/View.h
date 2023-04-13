#pragma once

#include "flan/Utility/Rect.h"

namespace flan {

/** Rectangular map from U to V. */
struct View
	{
	View( const Rect & _U, const Rect & _V )
		: U( _U )
		, V( _V )
		{}

	float wUToV( float w ) const { return w * V.w() / U.w(); }
	float hUToV( float h ) const { return h * V.h() / U.h(); }
	float wVToU( float w ) const { return w * U.w() / V.w(); }
	float hVToU( float h ) const { return h * U.h() / V.h(); }

	float xUToV( float x ) const { return wUToV( x - U.x1() ) + V.x1(); }
	float yUToV( float y ) const { return hUToV( y - U.y1() ) + V.y1(); }
	float xVToU( float x ) const { return wVToU( x - V.x1() ) + U.x1(); }
	float yVToU( float y ) const { return hVToU( y - V.y1() ) + U.y1(); }

	void forceAspectOnX( float aspect );
	void forceAspectOnY( float aspect );
	void growToAspect( float aspect );
	void shrinkToAspect( float aspect );

	Rect U, V;
	};

}