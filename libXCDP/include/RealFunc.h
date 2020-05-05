#pragma once

#include <functional>

#include "Utility.h"

namespace xcdp {

/*
 * Real function of one real variable
 */
struct RealFunc 
{
	template< typename Callable >
	RealFunc( Callable f_ ) : f( f_ ) {}
	RealFunc( float t0	  ) : f( [t0]( float ){ return t0; } ) {}
	RealFunc( double t0	  ) : f( [t0]( float ){ return t0; } ) {}
	RealFunc( int t0	  ) : f( [t0]( float ){ return float(t0); } ) {}
	RealFunc(			  ) : f( 0 ) {}

	float operator()( float t ) const { return f(t); }
	RealFunc operator()( const RealFunc g ) const { return [*this, g]( float t ){ return operator()( g(t) ); }; }
	float operator[]( float t ) const { return std::pow( 2.0f, f(t) ); }
    RealFunc logify() const { return [*this]( float t ){ return operator[](t); }; }

    RealFunc operator*( const RealFunc r ) const { return [*this, r](float t){ return operator()(t) * r(t); }; }
	RealFunc operator/( const RealFunc r ) const { return [*this, r](float t){ return operator()(t) / r(t); }; }
	RealFunc operator+( const RealFunc r ) const { return [*this, r](float t){ return operator()(t) + r(t); }; }
	RealFunc operator-( const RealFunc r ) const { return [*this, r](float t){ return operator()(t) - r(t); }; }
	RealFunc operator-() const { return [*this](float t){ return -operator()(t); }; }
	
	void graph( const std::string & filename,
		float left = -1, float right = 10, float bottom = -1, float top = 3.0, size_t resolution = 128 ) const;

	static RealFunc ADSR( float a, float d, float s, float r, float sLvl,
		float aExp = 1, float dExp = 1, float rExp = 1 );

	static RealFunc oscillate( RealFunc min, RealFunc max, RealFunc period = 2.0 * 3.14159265359, 
		RealFunc wave = []( float t ){ return sin(t); } );

	static RealFunc interpolatePoints( const std::vector< std::pair< float, float > > points, 
		Interpolator = Interpolators::linear );

private:
	std::function< float ( float ) > f;
};

struct Surface
{
	template< typename Callable >
	Surface( Callable f_ ) : f( f_ ) {}
	Surface( RealFunc f_ ) : f( [f_]( float t, float ){ return f_( t ); } ) {}
	Surface( float t0	 ) : f( [t0]( float, float ){ return t0; } ) {}
	Surface( int t0		 ) : f( [t0]( float, float ){ return float( t0 ); } ) {}
	float operator()( float time, float freq ) const { return f( time, freq ); }

	Surface operator*( const Surface r ) const { return [*this, r](float t, float f){ return operator()(t,f) * r(t,f); }; }
    Surface operator/( const Surface r ) const { return [*this, r](float t, float f){ return operator()(t,f) / r(t,f); }; }
    Surface operator+( const Surface r ) const { return [*this, r](float t, float f){ return operator()(t,f) + r(t,f); }; }
	Surface operator-( const Surface r ) const { return [*this, r](float t, float f){ return operator()(t,f) - r(t,f); }; }
	Surface operator-() const { return [*this](float t, float f){ return -operator()(t,f); }; }

private:
	std::function< float ( float, float ) > f;
};

} // End namespace xcdp
