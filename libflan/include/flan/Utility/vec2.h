#pragma once

#include <array>
#include <complex>

namespace flan {

/** Simple two element math vector. Arithmetic overloads are defined component-wise.
 */
struct vec2
	{
	float _x,_y;

	/** Default constructor, 0 initializes elements. */
	vec2() : _x( 0 ), _y( 0 ) {}

	/** Copy constructor */
	vec2( const vec2 & _v ) : _x( _v.x() ), _y( _v.y() ) {}

	/** Secondary copy constructor for easier syntax */
	vec2( float x, float y ) : _x( x ), _y( y ) {}

	/** Construct from std::complex */
	vec2( std::complex<float> c ) : _x( c.real() ), _y( c.imag() ) {}
	/** Convert to std::complex */
	operator std::complex<float>() const { return std::complex<float>( x(), y() ); }

	/** Construct from std::array */
	vec2( std::array<float,2> c ) : _x( c[0] ), _y( c[1] ) {}
	/** Convert to std::array */
	operator std::array<float,2>() const { return std::array<float,2>{ x(), y() }; }
 
	float & x() { return _x; }
	float & y() { return _y; }
	const float & x() const { return _x; }
	const float & y() const { return _y; }
	/** Identical to x() */
	float & t() { return _x; }
	/** Identical to y() */
	float & f() { return _y; }
	/** Identical to x() */
	const float & t() const { return _x; }
	/** Identical to y() */
	const float & f() const { return _y; }
	
	/** Scalar multiplication */
	vec2 operator*( const float s ) const  { return vec2{ x() * s, y() * s }; }
	/** Component-wise multiplication */
	vec2 operator*( const vec2 & r ) const  { return vec2{ x() * r.x(), y() * r.y() }; }
	/** Component-wise division */
	vec2 operator/( const vec2 & r ) const  { return vec2{ x() / r.x(), y() / r.y() }; }
	/** Component-wise addition */
	vec2 operator+( const vec2 & r ) const { return vec2{ x() + r.x(), y() + r.y() }; }
	/** Component-wise subtraction */
	vec2 operator-( const vec2 & r ) const { return vec2{ x() - r.x(), y() - r.y() }; }
	/** Component-wise negation */
	vec2 operator-() const { return vec2{ -x(), -y() }; }

	float mag() const { return std::sqrt( x() * x() + y() * y() ); }
	};

};